#!/usr/bin/env python3
"""Polish viral reference FASTAs from existing nf-core/rnaseq BAM outputs."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


VIRUS_SPECS = {
    "VEEV": {
        "datasets": ("mouse_veev", "rat_veev"),
        "reference_subdir": "VEEV",
    },
    "EEEV": {
        "datasets": ("mouse_eeev",),
        "reference_subdir": "EEEV",
    },
}


def default_viral_work_root(root_dir: Path) -> Path:
    preferred = root_dir / "viral_reference_work"
    legacy = root_dir / "viral_references"
    return preferred if preferred.exists() or not legacy.exists() else legacy


def parse_args() -> argparse.Namespace:
    root_dir = Path(__file__).resolve().parent.parent
    viral_work_root = default_viral_work_root(root_dir)
    parser = argparse.ArgumentParser(
        description=(
            "Reuse existing STAR BAM outputs to build guided viral consensus FASTAs "
            "for VEEV and/or EEEV."
        )
    )
    parser.add_argument(
        "virus",
        choices=("VEEV", "EEEV", "all"),
        help="Virus to polish, or 'all' to process both.",
    )
    parser.add_argument(
        "--results-root",
        default=str(root_dir / "results"),
        help="Root containing <dataset>/star_salmon/*.bam outputs (default: repo results/).",
    )
    parser.add_argument(
        "--reference-root",
        default=str(root_dir / "references"),
        help="Reference root containing VEEV/virus.fa and EEEV/virus.fa (default: repo references/).",
    )
    parser.add_argument(
        "--output-root",
        default=str(viral_work_root / "polish"),
        help="Output root for polished FASTAs, VCFs, and summaries.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of top-ranked BAMs to use per virus (default: 5).",
    )
    parser.add_argument(
        "--min-mapped-reads",
        type=int,
        default=0,
        help="Minimum viral mapped reads required to keep a BAM candidate (default: 0).",
    )
    parser.add_argument(
        "--mask-depth",
        type=int,
        default=10,
        help="Mask consensus positions below this total depth with N (default: 10).",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=1_000_000,
        help="bcftools mpileup max depth cap (default: 1000000).",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=4,
        help="Concurrent worker count for per-BAM ranking and per-sample consensus generation (default: 4).",
    )
    parser.add_argument(
        "--samtools-bin",
        default=os.environ.get("VIRAL_SAMTOOLS_BIN", "samtools"),
        help="samtools executable to use (default: env VIRAL_SAMTOOLS_BIN or samtools).",
    )
    parser.add_argument(
        "--bcftools-bin",
        default=os.environ.get("VIRAL_BCFTOOLS_BIN", "bcftools"),
        help="bcftools executable to use (default: env VIRAL_BCFTOOLS_BIN or bcftools).",
    )
    parser.add_argument(
        "--install",
        action="store_true",
        help="Overwrite references/<virus>/virus.fa with the pooled polished FASTA after backup.",
    )
    return parser.parse_args()


def require_command(path_or_name: str) -> str:
    resolved = shutil.which(path_or_name)
    if resolved:
        return resolved
    candidate = Path(path_or_name).expanduser()
    if candidate.is_file():
        return str(candidate)
    raise FileNotFoundError(f"Required executable not found: {path_or_name}")


def read_single_fasta_header(ref_fasta: Path) -> tuple[str, int]:
    seqid = None
    length = 0
    with ref_fasta.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seqid is not None:
                    raise ValueError(f"Expected a single FASTA record in {ref_fasta}")
                seqid = line[1:].split()[0]
                continue
            length += len(line)
    if seqid is None:
        raise ValueError(f"No FASTA header found in {ref_fasta}")
    return seqid, length


def ensure_fasta_index(samtools_bin: str, ref_fasta: Path) -> None:
    fai_path = ref_fasta.with_suffix(ref_fasta.suffix + ".fai")
    if fai_path.is_file():
        return
    subprocess.run([samtools_bin, "faidx", str(ref_fasta)], check=True, capture_output=True, text=True)


def discover_bams(results_root: Path, datasets: tuple[str, ...]) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    for dataset in datasets:
        bam_dir = results_root / dataset / "star_salmon"
        if not bam_dir.is_dir():
            raise FileNotFoundError(f"STAR BAM directory not found: {bam_dir}")

        for bam_path in sorted(bam_dir.glob("*.bam")):
            if bam_path.name.endswith(".bam.bai"):
                continue
            sample_id = bam_path.name
            if sample_id.endswith(".markdup.sorted.bam"):
                sample_id = sample_id[: -len(".markdup.sorted.bam")]
            else:
                sample_id = bam_path.stem
            candidates.append(
                {
                    "dataset": dataset,
                    "sample_id": sample_id,
                    "bam_path": str(bam_path.resolve()),
                }
            )
    if not candidates:
        raise RuntimeError(f"No BAM candidates found under {results_root}")
    return candidates


def parse_idxstats_for_contig(samtools_bin: str, bam_path: Path, contig: str) -> tuple[int, int, int]:
    result = subprocess.run(
        [samtools_bin, "idxstats", str(bam_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    for line in result.stdout.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 4:
            continue
        seqid, length, mapped, unmapped = fields
        if seqid == contig:
            return int(length), int(mapped), int(unmapped)
    raise ValueError(f"Contig {contig} was not found in {bam_path}")


def summarize_depth(
    samtools_bin: str,
    bam_paths: list[Path],
    contig: str,
    genome_length: int,
    mask_depth: int,
) -> tuple[float, float]:
    cmd = [samtools_bin, "depth", "-aa", "-r", contig, *map(str, bam_paths)]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)

    total_depth = 0
    covered_bases = 0
    observed_rows = 0
    for line in result.stdout.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        depth_total = sum(int(value) for value in fields[2:])
        total_depth += depth_total
        if depth_total >= mask_depth:
            covered_bases += 1
        observed_rows += 1

    if observed_rows == 0:
        return 0.0, 0.0

    mean_depth = total_depth / observed_rows
    breadth = covered_bases / genome_length if genome_length > 0 else 0.0
    return mean_depth, breadth


def write_low_depth_mask(
    samtools_bin: str,
    bam_paths: list[Path],
    contig: str,
    mask_depth: int,
    bed_path: Path,
) -> tuple[int, int]:
    cmd = [samtools_bin, "depth", "-aa", "-r", contig, *map(str, bam_paths)]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)

    masked_bases = 0
    masked_regions = 0
    start_zero = None
    prev_pos = None

    bed_path.parent.mkdir(parents=True, exist_ok=True)
    with bed_path.open("w", encoding="utf-8") as handle:
        for line in result.stdout.splitlines():
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            pos = int(fields[1])
            depth_total = sum(int(value) for value in fields[2:])
            if depth_total < mask_depth:
                masked_bases += 1
                if start_zero is None:
                    start_zero = pos - 1
                prev_pos = pos
                continue

            if start_zero is not None and prev_pos is not None:
                handle.write(f"{contig}\t{start_zero}\t{prev_pos}\n")
                masked_regions += 1
                start_zero = None
                prev_pos = None

        if start_zero is not None and prev_pos is not None:
            handle.write(f"{contig}\t{start_zero}\t{prev_pos}\n")
            masked_regions += 1

    return masked_bases, masked_regions


def build_variant_calls(
    bcftools_bin: str,
    ref_fasta: Path,
    contig: str,
    bam_paths: list[Path],
    max_depth: int,
    vcf_path: Path,
) -> None:
    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    mpileup = subprocess.Popen(
        [
            bcftools_bin,
            "mpileup",
            "-Ou",
            "-f",
            str(ref_fasta),
            "-r",
            contig,
            "-d",
            str(max_depth),
            *map(str, bam_paths),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
    )
    caller = subprocess.Popen(
        [
            bcftools_bin,
            "call",
            "--ploidy",
            "1",
            "-mv",
            "-Oz",
            "-o",
            str(vcf_path),
        ],
        stdin=mpileup.stdout,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=False,
    )
    assert mpileup.stdout is not None
    mpileup.stdout.close()
    caller_stderr = caller.communicate()[1]
    mpileup_stderr = mpileup.communicate()[1]

    if mpileup.returncode != 0:
        raise RuntimeError(mpileup_stderr.decode("utf-8", errors="replace").strip() or "bcftools mpileup failed")
    if caller.returncode != 0:
        raise RuntimeError(caller_stderr.decode("utf-8", errors="replace").strip() or "bcftools call failed")

    subprocess.run([bcftools_bin, "index", "-f", str(vcf_path)], check=True, capture_output=True, text=True)


def count_vcf_variants(bcftools_bin: str, vcf_path: Path) -> int:
    result = subprocess.run(
        [bcftools_bin, "view", "-H", str(vcf_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return sum(1 for line in result.stdout.splitlines() if line.strip())


def write_consensus_fasta(
    bcftools_bin: str,
    ref_fasta: Path,
    vcf_path: Path,
    mask_bed: Path,
    output_fasta: Path,
) -> None:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [
            bcftools_bin,
            "consensus",
            "-f",
            str(ref_fasta),
            "--mask",
            str(mask_bed),
            str(vcf_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    output_fasta.write_text(result.stdout, encoding="utf-8")


def write_summary(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def collect_candidate_metrics(
    samtools_bin: str,
    bam_path: Path,
    dataset: str,
    sample_id: str,
    virus_id: str,
    contig: str,
    genome_length: int,
    mask_depth: int,
) -> dict[str, object]:
    _length, mapped, unmapped = parse_idxstats_for_contig(samtools_bin, bam_path, contig)
    mean_depth, breadth = summarize_depth(
        samtools_bin,
        [bam_path],
        contig,
        genome_length,
        mask_depth,
    )
    return {
        "virus_id": virus_id,
        "dataset": dataset,
        "sample_id": sample_id,
        "bam_path": str(bam_path),
        "viral_contig": contig,
        "viral_mapped_reads": mapped,
        "viral_unmapped_reads": unmapped,
        "mean_depth": f"{mean_depth:.4f}",
        "breadth_ge_mask_depth": f"{breadth:.6f}",
        "selected": "no",
        "variant_count": "",
        "masked_bases": "",
        "masked_regions": "",
        "consensus_fasta": "",
        "vcf_path": "",
        "mask_bed": "",
    }


def build_sample_consensus(
    args: argparse.Namespace,
    ref_fasta: Path,
    contig: str,
    samples_dir: Path,
    row: dict[str, object],
) -> dict[str, object]:
    row = dict(row)
    row["selected"] = "yes"
    sample_prefix = samples_dir / str(row["sample_id"])
    sample_bam = Path(str(row["bam_path"]))
    sample_vcf = sample_prefix.with_suffix(".vcf.gz")
    sample_mask = sample_prefix.with_suffix(".lowcov.bed")
    sample_consensus = sample_prefix.with_suffix(".consensus.fa")

    build_variant_calls(
        args.bcftools_bin,
        ref_fasta,
        contig,
        [sample_bam],
        args.max_depth,
        sample_vcf,
    )
    masked_bases, masked_regions = write_low_depth_mask(
        args.samtools_bin,
        [sample_bam],
        contig,
        args.mask_depth,
        sample_mask,
    )
    write_consensus_fasta(
        args.bcftools_bin,
        ref_fasta,
        sample_vcf,
        sample_mask,
        sample_consensus,
    )

    row["variant_count"] = count_vcf_variants(args.bcftools_bin, sample_vcf)
    row["masked_bases"] = masked_bases
    row["masked_regions"] = masked_regions
    row["consensus_fasta"] = str(sample_consensus)
    row["vcf_path"] = str(sample_vcf)
    row["mask_bed"] = str(sample_mask)
    return row


def polish_one_virus(args: argparse.Namespace, virus_id: str) -> None:
    spec = VIRUS_SPECS[virus_id]
    results_root = Path(args.results_root).expanduser().resolve()
    reference_root = Path(args.reference_root).expanduser().resolve()
    output_root = Path(args.output_root).expanduser().resolve()

    ref_dir = reference_root / spec["reference_subdir"]
    ref_fasta = ref_dir / "virus.fa"
    if not ref_fasta.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")

    contig, genome_length = read_single_fasta_header(ref_fasta)
    ensure_fasta_index(args.samtools_bin, ref_fasta)

    bam_rows = discover_bams(results_root, spec["datasets"])
    ranked_rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = [
            pool.submit(
                collect_candidate_metrics,
                args.samtools_bin,
                Path(str(row["bam_path"])),
                str(row["dataset"]),
                str(row["sample_id"]),
                virus_id,
                contig,
                genome_length,
                args.mask_depth,
            )
            for row in bam_rows
        ]
        for future in as_completed(futures):
            ranked_rows.append(future.result())

    ranked_rows.sort(
        key=lambda row: (
            int(row["viral_mapped_reads"]),
            float(row["mean_depth"]),
            str(row["sample_id"]),
        ),
        reverse=True,
    )

    selected_rows = [
        row
        for row in ranked_rows
        if int(row["viral_mapped_reads"]) > 0 and int(row["viral_mapped_reads"]) >= args.min_mapped_reads
    ][: args.top_n]
    if not selected_rows:
        raise RuntimeError(
            f"No BAMs met the selection criteria for {virus_id}. "
            f"Try lowering --min-mapped-reads or checking the results root."
        )

    virus_output_dir = output_root / virus_id
    samples_dir = virus_output_dir / "samples"
    samples_dir.mkdir(parents=True, exist_ok=True)

    selected_bams = [Path(str(row["bam_path"])) for row in selected_rows]
    selected_by_sample = {}
    with ThreadPoolExecutor(max_workers=max(1, min(args.jobs, len(selected_rows)))) as pool:
        futures = [
            pool.submit(
                build_sample_consensus,
                args,
                ref_fasta,
                contig,
                samples_dir,
                row,
            )
            for row in selected_rows
        ]
        for future in as_completed(futures):
            finished = future.result()
            selected_by_sample[str(finished["sample_id"])] = finished

    ranked_rows = [
        selected_by_sample.get(str(row["sample_id"]), row)
        for row in ranked_rows
    ]
    selected_rows = [selected_by_sample[str(row["sample_id"])] for row in selected_rows]

    pooled_vcf = virus_output_dir / "pooled.vcf.gz"
    pooled_mask = virus_output_dir / "pooled.lowcov.bed"
    pooled_fasta = virus_output_dir / "virus.fa"
    pooled_selected = virus_output_dir / "selected_bams.txt"

    build_variant_calls(
        args.bcftools_bin,
        ref_fasta,
        contig,
        selected_bams,
        args.max_depth,
        pooled_vcf,
    )
    pooled_masked_bases, pooled_masked_regions = write_low_depth_mask(
        args.samtools_bin,
        selected_bams,
        contig,
        args.mask_depth,
        pooled_mask,
    )
    write_consensus_fasta(
        args.bcftools_bin,
        ref_fasta,
        pooled_vcf,
        pooled_mask,
        pooled_fasta,
    )

    pooled_mean_depth, pooled_breadth = summarize_depth(
        args.samtools_bin,
        selected_bams,
        contig,
        genome_length,
        args.mask_depth,
    )
    pooled_summary = [
        {
            "virus_id": virus_id,
            "viral_contig": contig,
            "genome_length": genome_length,
            "selected_bam_count": len(selected_bams),
            "selected_samples": ",".join(str(row["sample_id"]) for row in selected_rows),
            "pooled_variant_count": count_vcf_variants(args.bcftools_bin, pooled_vcf),
            "pooled_masked_bases": pooled_masked_bases,
            "pooled_masked_regions": pooled_masked_regions,
            "pooled_mean_depth": f"{pooled_mean_depth:.4f}",
            "pooled_breadth_ge_mask_depth": f"{pooled_breadth:.6f}",
            "reference_fasta": str(ref_fasta),
            "pooled_vcf": str(pooled_vcf),
            "pooled_mask_bed": str(pooled_mask),
            "pooled_consensus_fasta": str(pooled_fasta),
            "installed_reference_fasta": str(ref_fasta if args.install else ""),
        }
    ]

    with pooled_selected.open("w", encoding="utf-8") as handle:
        for bam_path in selected_bams:
            handle.write(f"{bam_path}\n")

    write_summary(
        virus_output_dir / "sample_summary.tsv",
        ranked_rows,
        [
            "virus_id",
            "dataset",
            "sample_id",
            "bam_path",
            "viral_contig",
            "viral_mapped_reads",
            "viral_unmapped_reads",
            "mean_depth",
            "breadth_ge_mask_depth",
            "selected",
            "variant_count",
            "masked_bases",
            "masked_regions",
            "consensus_fasta",
            "vcf_path",
            "mask_bed",
        ],
    )
    write_summary(
        virus_output_dir / "pooled_summary.tsv",
        pooled_summary,
        [
            "virus_id",
            "viral_contig",
            "genome_length",
            "selected_bam_count",
            "selected_samples",
            "pooled_variant_count",
            "pooled_masked_bases",
            "pooled_masked_regions",
            "pooled_mean_depth",
            "pooled_breadth_ge_mask_depth",
            "reference_fasta",
            "pooled_vcf",
            "pooled_mask_bed",
            "pooled_consensus_fasta",
            "installed_reference_fasta",
        ],
    )

    if args.install:
        backup_path = virus_output_dir / "virus.pre_polish.fa"
        shutil.copyfile(ref_fasta, backup_path)
        shutil.copyfile(pooled_fasta, ref_fasta)
        subprocess.run([args.samtools_bin, "faidx", str(ref_fasta)], check=True, capture_output=True, text=True)

    print(
        f"{virus_id}: polished consensus written to {pooled_fasta} "
        f"from {len(selected_bams)} BAM(s) on contig {contig}",
        file=sys.stderr,
    )


def main() -> int:
    args = parse_args()
    if args.top_n <= 0:
        raise ValueError("--top-n must be a positive integer")
    if args.jobs <= 0:
        raise ValueError("--jobs must be a positive integer")
    if args.min_mapped_reads < 0:
        raise ValueError("--min-mapped-reads cannot be negative")
    if args.mask_depth <= 0:
        raise ValueError("--mask-depth must be a positive integer")
    if args.max_depth <= 0:
        raise ValueError("--max-depth must be a positive integer")
    args.samtools_bin = require_command(args.samtools_bin)
    args.bcftools_bin = require_command(args.bcftools_bin)

    viruses = ("VEEV", "EEEV") if args.virus == "all" else (args.virus,)
    for virus_id in viruses:
        polish_one_virus(args, virus_id)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
