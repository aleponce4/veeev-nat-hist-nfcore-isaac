#!/usr/bin/env python3
"""Empirically derive alphavirus 26S transcript starts from RNA-seq data."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


PHASE1_PAIRS_DEFAULT = 100_000
PHASE2_PAIRS_DEFAULT = 1_000_000
TOP_N_DEFAULT = 5
SEARCH_UPSTREAM_DEFAULT = 400
SEARCH_DOWNSTREAM_DEFAULT = 200
JUMP_WINDOW_DEFAULT = 10
JUMP_MULTIPLIER_DEFAULT = 5.0
JOBS_DEFAULT = 4


VIRUS_SPECS = {
    "VEEV": {
        "display_name": "VEEV_TrD",
        "source_fasta": "VEEV/VEEV_TrD.fa",
        "source_cds_gtf": "VEEV/VEEV_TrD (CDS).gtf",
        "candidate_datasets": ("mouse_veev", "rat_veev"),
        "gene_id": "VEEV",
        "tx_49s": "VEEV_49S",
        "tx_26s": "VEEV_26S",
        "gene_name": "VEEV",
    },
    "EEEV": {
        "display_name": "EEEV_FL93",
        "source_fasta": "EEEV/EEEV_FL93.fa",
        "source_cds_gtf": "EEEV/EEEV_FL93 (CDS).gtf",
        "candidate_datasets": ("mouse_eeev",),
        "gene_id": "EEEV",
        "tx_49s": "EEEV_49S",
        "tx_26s": "EEEV_26S",
        "gene_name": "EEEV",
    },
}


def parse_args() -> argparse.Namespace:
    root_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description=(
            "Use a two-tier viral-only pilot mapping workflow to derive empirical "
            "26S transcript starts and write curated viral FASTA/GTF files."
        )
    )
    parser.add_argument(
        "virus",
        choices=("VEEV", "EEEV", "all"),
        help="Virus to curate, or 'all' to build both VEEV and EEEV references.",
    )
    parser.add_argument(
        "--raw-root",
        default=str(root_dir / "viral_references" / "raw"),
        help="Tracked raw viral reference root (default: repo viral_references/raw)",
    )
    parser.add_argument(
        "--input-root",
        default=str(root_dir / "inputs"),
        help="Input dataset root containing mouse_veev, mouse_eeev, and rat_veev directories.",
    )
    parser.add_argument(
        "--output-root",
        default=str(root_dir / "references"),
        help="Reference output root containing VEEV/ and EEEV/ directories.",
    )
    parser.add_argument(
        "--results-dir",
        default=str(root_dir / "viral_references" / "curation"),
        help="Directory for tracked empirical TSS metadata tables.",
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="Optional staging manifest CSV. If provided, sample metadata will be read from it.",
    )
    parser.add_argument(
        "--phase1-pairs",
        type=int,
        default=PHASE1_PAIRS_DEFAULT,
        help=f"Read pairs to map for the all-sample ranking phase (default: {PHASE1_PAIRS_DEFAULT}).",
    )
    parser.add_argument(
        "--phase2-pairs",
        type=int,
        default=PHASE2_PAIRS_DEFAULT,
        help=f"Read pairs to map for the top-sample deep phase (default: {PHASE2_PAIRS_DEFAULT}).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=TOP_N_DEFAULT,
        help=f"Number of top-ranked samples to deep-map per virus (default: {TOP_N_DEFAULT}).",
    )
    parser.add_argument(
        "--search-upstream",
        type=int,
        default=SEARCH_UPSTREAM_DEFAULT,
        help=f"Bases upstream of the first structural CDS to scan for the 26S TSS (default: {SEARCH_UPSTREAM_DEFAULT}).",
    )
    parser.add_argument(
        "--search-downstream",
        type=int,
        default=SEARCH_DOWNSTREAM_DEFAULT,
        help=f"Bases downstream of the first structural CDS to scan for the 26S TSS (default: {SEARCH_DOWNSTREAM_DEFAULT}).",
    )
    parser.add_argument(
        "--jump-window",
        type=int,
        default=JUMP_WINDOW_DEFAULT,
        help=f"Window size for the previous-depth average used in slope detection (default: {JUMP_WINDOW_DEFAULT}).",
    )
    parser.add_argument(
        "--jump-multiplier",
        type=float,
        default=JUMP_MULTIPLIER_DEFAULT,
        help=f"Required fold jump over the previous-depth average (default: {JUMP_MULTIPLIER_DEFAULT}).",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=JOBS_DEFAULT,
        help=f"Number of concurrent pilot jobs (default: {JOBS_DEFAULT}).",
    )
    parser.add_argument(
        "--mapper-bin",
        default=os.environ.get("VIRAL_MAPPER_BIN", "minimap2"),
        help="Mapper executable to use for viral-only pilot alignments (default: env VIRAL_MAPPER_BIN or minimap2).",
    )
    parser.add_argument(
        "--samtools-bin",
        default=os.environ.get("VIRAL_SAMTOOLS_BIN", "samtools"),
        help="samtools executable to use (default: env VIRAL_SAMTOOLS_BIN or samtools).",
    )
    parser.add_argument(
        "--gene-id",
        default=None,
        help="Optional override for the curated viral gene_id. Defaults to the virus name (VEEV or EEEV).",
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


def parse_manifest(manifest_path: Path) -> dict[str, list[dict[str, str]]]:
    rows_by_dataset: dict[str, list[dict[str, str]]] = {}
    with manifest_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows_by_dataset.setdefault(row["dataset"], []).append(row)
    return rows_by_dataset


def collect_samples(input_root: Path, datasets: tuple[str, ...], manifest_path: Path | None) -> list[dict[str, str]]:
    rows_by_dataset = parse_manifest(manifest_path) if manifest_path else {}
    samples: list[dict[str, str]] = []

    for dataset in datasets:
        if manifest_path:
            dataset_rows = rows_by_dataset.get(dataset, [])
            for row in dataset_rows:
                samples.append(
                    {
                        "dataset": dataset,
                        "sample_id": row["sample_id"],
                        "fastq_1": row["staged_fastq_1"] or row["source_fastq_1"],
                        "fastq_2": row["staged_fastq_2"] or row["source_fastq_2"],
                    }
                )
            continue

        dataset_dir = input_root / dataset
        if not dataset_dir.is_dir():
            raise FileNotFoundError(f"Input dataset directory not found: {dataset_dir}")

        for r1_path in sorted(dataset_dir.glob("*_R1_001.fastq.gz")):
            sample_id = r1_path.name.removesuffix("_R1_001.fastq.gz")
            r2_path = dataset_dir / f"{sample_id}_R2_001.fastq.gz"
            if not r2_path.exists():
                raise FileNotFoundError(f"Missing R2 FASTQ for sample {sample_id}: {r2_path}")
            samples.append(
                {
                    "dataset": dataset,
                    "sample_id": sample_id,
                    "fastq_1": str(r1_path.resolve()),
                    "fastq_2": str(r2_path.resolve()),
                }
            )

    if not samples:
        raise RuntimeError(f"No candidate FASTQ pairs found for datasets: {', '.join(datasets)}")
    return samples


def subset_fastq_pair(
    r1_path: Path,
    r2_path: Path,
    r1_subset: Path,
    r2_subset: Path,
    pair_limit: int,
) -> int:
    sampled_pairs = 0
    with gzip.open(r1_path, "rt", encoding="utf-8") as r1_handle, gzip.open(
        r2_path, "rt", encoding="utf-8"
    ) as r2_handle, r1_subset.open("w", encoding="utf-8") as r1_out, r2_subset.open(
        "w", encoding="utf-8"
    ) as r2_out:
        while sampled_pairs < pair_limit:
            r1_chunk = [r1_handle.readline() for _ in range(4)]
            r2_chunk = [r2_handle.readline() for _ in range(4)]
            if not all(r1_chunk) or not all(r2_chunk):
                break
            r1_out.writelines(r1_chunk)
            r2_out.writelines(r2_chunk)
            sampled_pairs += 1
    if sampled_pairs == 0:
        raise RuntimeError(f"No reads were sampled from {r1_path} / {r2_path}")
    return sampled_pairs


def run_mapped_primary_count(
    mapper_bin: str,
    samtools_bin: str,
    ref_fasta: Path,
    r1_subset: Path,
    r2_subset: Path,
    stderr_path: Path,
) -> int:
    with stderr_path.open("w", encoding="utf-8") as stderr_handle:
        mapper = subprocess.Popen(
            [mapper_bin, "-ax", "sr", str(ref_fasta), str(r1_subset), str(r2_subset)],
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            text=False,
        )
        counter = subprocess.Popen(
            [samtools_bin, "view", "-c", "-F", "2308", "-"],
            stdin=mapper.stdout,
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            text=True,
        )
        assert mapper.stdout is not None
        mapper.stdout.close()
        stdout, _ = counter.communicate()
        mapper_return = mapper.wait()

    if mapper_return != 0:
        raise RuntimeError(f"Mapper failed for subset FASTQs; see {stderr_path}")
    if counter.returncode != 0:
        raise RuntimeError(f"samtools view failed for subset FASTQs; see {stderr_path}")
    return int(stdout.strip() or "0")


def run_sorted_bam(
    mapper_bin: str,
    samtools_bin: str,
    ref_fasta: Path,
    r1_subset: Path,
    r2_subset: Path,
    bam_path: Path,
    stderr_path: Path,
) -> None:
    with stderr_path.open("w", encoding="utf-8") as stderr_handle:
        mapper = subprocess.Popen(
            [mapper_bin, "-ax", "sr", str(ref_fasta), str(r1_subset), str(r2_subset)],
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            text=False,
        )
        sorter = subprocess.Popen(
            [samtools_bin, "sort", "-o", str(bam_path), "-"],
            stdin=mapper.stdout,
            stdout=subprocess.DEVNULL,
            stderr=stderr_handle,
            text=False,
        )
        assert mapper.stdout is not None
        mapper.stdout.close()
        sorter.wait()
        mapper_return = mapper.wait()

    if mapper_return != 0 or sorter.returncode != 0:
        raise RuntimeError(f"Pilot BAM creation failed; see {stderr_path}")

    subprocess.run([samtools_bin, "index", str(bam_path)], check=True, capture_output=True, text=True)


def read_fasta_header_and_length(fasta_path: Path) -> tuple[str, int]:
    seqid = None
    length = 0
    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seqid is not None:
                    raise ValueError(f"Expected a single FASTA record in {fasta_path}")
                seqid = line[1:].split()[0]
                continue
            length += len(line)
    if seqid is None:
        raise ValueError(f"No FASTA header found in {fasta_path}")
    return seqid, length


def detect_structural_start(cds_gtf_path: Path) -> int:
    structural_starts = []
    with cds_gtf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            seqid, _source, feature, start, _end, _score, _strand, _frame, attrs = line.rstrip("\n").split("\t")
            if feature != "CDS":
                continue
            gene_match = attrs.split('gene_id "', 1)
            if len(gene_match) != 2:
                continue
            gene_name = gene_match[1].split('"', 1)[0]
            upper_gene = gene_name.upper()
            if upper_gene.startswith("NSP") or upper_gene.startswith("CDS_"):
                continue
            structural_starts.append(int(start))
    if not structural_starts:
        raise ValueError(f"Could not identify structural CDS starts in {cds_gtf_path}")
    return min(structural_starts)


def detect_tss_from_depth(
    samtools_bin: str,
    bam_path: Path,
    seqid: str,
    search_start: int,
    search_end: int,
    jump_window: int,
    jump_multiplier: float,
) -> tuple[int, list[tuple[int, int]]]:
    result = subprocess.run(
        [samtools_bin, "depth", "-aa", "-r", f"{seqid}:{search_start}-{search_end}", str(bam_path)],
        check=True,
        capture_output=True,
        text=True,
    )

    depth_rows: list[tuple[int, int]] = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        _seqid, pos, depth = line.split("\t")
        depth_rows.append((int(pos), int(depth)))

    if len(depth_rows) <= jump_window:
        raise ValueError(f"Depth window is too short for TSS detection in {bam_path}")

    best_ratio = -math.inf
    best_pos = None
    for idx in range(jump_window, len(depth_rows)):
        current_pos, current_depth = depth_rows[idx]
        prev_depths = [depth for _pos, depth in depth_rows[idx - jump_window : idx]]
        prev_avg = sum(prev_depths) / jump_window
        if prev_avg <= 0:
            continue
        ratio = current_depth / prev_avg
        if ratio > best_ratio:
            best_ratio = ratio
            best_pos = current_pos
        if ratio >= jump_multiplier:
            return current_pos, depth_rows

    if best_pos is None:
        raise ValueError(f"No usable coverage jump was found in {bam_path}")
    return best_pos, depth_rows


def build_curated_gtf(
    seqid: str,
    genome_length: int,
    consensus_tss: int,
    gene_id: str,
    gene_name: str,
    tx_49s: str,
    tx_26s: str,
) -> str:
    if not (1 <= consensus_tss <= genome_length):
        raise ValueError(f"Consensus TSS {consensus_tss} is outside the genome length {genome_length}")

    attrs_gene = f'gene_id "{gene_id}"; gene_name "{gene_name}";'
    attrs_49s = (
        f'gene_id "{gene_id}"; gene_name "{gene_name}"; '
        f'transcript_id "{tx_49s}"; transcript_name "{tx_49s}";'
    )
    attrs_26s = (
        f'gene_id "{gene_id}"; gene_name "{gene_name}"; '
        f'transcript_id "{tx_26s}"; transcript_name "{tx_26s}";'
    )

    rows = [
        (seqid, "curated", "gene", 1, genome_length, ".", "+", ".", attrs_gene),
        (seqid, "curated", "transcript", 1, genome_length, ".", "+", ".", attrs_49s),
        (seqid, "curated", "exon", 1, genome_length, ".", "+", ".", attrs_49s),
        (seqid, "curated", "transcript", consensus_tss, genome_length, ".", "+", ".", attrs_26s),
        (seqid, "curated", "exon", consensus_tss, genome_length, ".", "+", ".", attrs_26s),
    ]
    return "\n".join("\t".join(str(field) for field in row) for row in rows) + "\n"


def analyze_phase1_sample(
    sample: dict[str, str],
    mapper_bin: str,
    samtools_bin: str,
    ref_fasta: Path,
    phase1_pairs: int,
    temp_root: Path,
) -> dict[str, object]:
    sample_id = sample["sample_id"]
    dataset = sample["dataset"]
    work_dir = temp_root / f"phase1_{dataset}_{sample_id}"
    work_dir.mkdir(parents=True, exist_ok=True)
    r1_subset = work_dir / "subset_R1.fastq"
    r2_subset = work_dir / "subset_R2.fastq"
    sampled_pairs = subset_fastq_pair(Path(sample["fastq_1"]), Path(sample["fastq_2"]), r1_subset, r2_subset, phase1_pairs)
    mapped_primary_reads = run_mapped_primary_count(
        mapper_bin,
        samtools_bin,
        ref_fasta,
        r1_subset,
        r2_subset,
        work_dir / "phase1.stderr.log",
    )
    total_reads_sampled = sampled_pairs * 2
    viral_fraction = mapped_primary_reads / total_reads_sampled
    return {
        "sample_id": sample_id,
        "dataset": dataset,
        "fastq_1": sample["fastq_1"],
        "fastq_2": sample["fastq_2"],
        "phase1_pairs": sampled_pairs,
        "phase1_total_reads": total_reads_sampled,
        "phase1_mapped_primary_reads": mapped_primary_reads,
        "phase1_viral_fraction": viral_fraction,
    }


def analyze_phase2_sample(
    sample: dict[str, object],
    mapper_bin: str,
    samtools_bin: str,
    ref_fasta: Path,
    seqid: str,
    search_start: int,
    search_end: int,
    phase2_pairs: int,
    jump_window: int,
    jump_multiplier: float,
    temp_root: Path,
) -> tuple[dict[str, object], list[tuple[int, int]]]:
    sample_id = str(sample["sample_id"])
    dataset = str(sample["dataset"])
    work_dir = temp_root / f"phase2_{dataset}_{sample_id}"
    work_dir.mkdir(parents=True, exist_ok=True)
    r1_subset = work_dir / "subset_R1.fastq"
    r2_subset = work_dir / "subset_R2.fastq"
    sampled_pairs = subset_fastq_pair(
        Path(str(sample["fastq_1"])),
        Path(str(sample["fastq_2"])),
        r1_subset,
        r2_subset,
        phase2_pairs,
    )
    bam_path = work_dir / "pilot_virus.bam"
    run_sorted_bam(
        mapper_bin,
        samtools_bin,
        ref_fasta,
        r1_subset,
        r2_subset,
        bam_path,
        work_dir / "phase2.stderr.log",
    )
    tss_coord, depth_rows = detect_tss_from_depth(
        samtools_bin,
        bam_path,
        seqid,
        search_start,
        search_end,
        jump_window,
        jump_multiplier,
    )
    sample = dict(sample)
    sample["phase2_pairs"] = sampled_pairs
    sample["phase2_tss_coord"] = tss_coord
    return sample, depth_rows


def write_results_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "virus_id",
        "dataset",
        "sample_id",
        "phase1_pairs",
        "phase1_total_reads",
        "phase1_mapped_primary_reads",
        "phase1_viral_fraction",
        "selected_top_n",
        "phase2_status",
        "phase2_pairs",
        "phase2_tss_coord",
        "phase2_error",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_consensus_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "virus_id",
        "source_fasta",
        "source_cds_gtf",
        "candidate_datasets",
        "selected_samples",
        "seqid",
        "genome_length",
        "structural_start",
        "search_start",
        "search_end",
        "phase1_pairs",
        "phase2_pairs",
        "top_n",
        "jump_window",
        "jump_multiplier",
        "consensus_tss",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def curate_one_virus(args: argparse.Namespace, virus_id: str) -> tuple[list[dict[str, object]], dict[str, object]]:
    spec = VIRUS_SPECS[virus_id]
    raw_root = Path(args.raw_root)
    input_root = Path(args.input_root)
    output_root = Path(args.output_root)

    ref_fasta = raw_root / spec["source_fasta"]
    source_cds_gtf = raw_root / spec["source_cds_gtf"]
    if not ref_fasta.is_file():
        raise FileNotFoundError(f"Raw viral FASTA not found: {ref_fasta}")
    if not source_cds_gtf.is_file():
        raise FileNotFoundError(f"Raw viral CDS GTF not found: {source_cds_gtf}")

    seqid, genome_length = read_fasta_header_and_length(ref_fasta)
    structural_start = detect_structural_start(source_cds_gtf)
    search_start = max(1, structural_start - args.search_upstream)
    search_end = min(genome_length, structural_start + args.search_downstream)

    manifest_path = Path(args.manifest).expanduser().resolve() if args.manifest else None
    samples = collect_samples(input_root, spec["candidate_datasets"], manifest_path)

    with tempfile.TemporaryDirectory(prefix=f"{virus_id.lower()}_tss_") as temp_dir:
        temp_root = Path(temp_dir)
        phase1_rows: list[dict[str, object]] = []
        print(
            f"[{virus_id}] Phase 1: ranking {len(samples)} samples with {args.phase1_pairs} read pairs each",
            file=sys.stderr,
            flush=True,
        )
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = [
                pool.submit(
                    analyze_phase1_sample,
                    sample,
                    args.mapper_bin,
                    args.samtools_bin,
                    ref_fasta,
                    args.phase1_pairs,
                    temp_root,
                )
                for sample in samples
            ]
            for future in as_completed(futures):
                phase1_rows.append(future.result())

        phase1_rows.sort(
            key=lambda row: (
                float(row["phase1_viral_fraction"]),
                int(row["phase1_mapped_primary_reads"]),
                row["sample_id"],
            ),
            reverse=True,
        )

        target_successes = min(args.top_n, len(phase1_rows))
        print(
            f"[{virus_id}] Phase 2: deep-mapping up to {target_successes} top samples with {args.phase2_pairs} read pairs each",
            file=sys.stderr,
            flush=True,
        )
        phase2_rows: list[dict[str, object]] = []
        phase2_failures: dict[str, str] = {}
        for row in phase1_rows:
            if len(phase2_rows) >= target_successes:
                break
            sample_label = f'{row["dataset"]}/{row["sample_id"]}'
            try:
                phase2_row, _depth_rows = analyze_phase2_sample(
                    row,
                    args.mapper_bin,
                    args.samtools_bin,
                    ref_fasta,
                    seqid,
                    search_start,
                    search_end,
                    args.phase2_pairs,
                    args.jump_window,
                    args.jump_multiplier,
                    temp_root,
                )
                phase2_rows.append(phase2_row)
                print(
                    f'[{virus_id}] Phase 2 success: {sample_label} -> TSS {phase2_row["phase2_tss_coord"]}',
                    file=sys.stderr,
                    flush=True,
                )
            except Exception as exc:  # noqa: BLE001
                phase2_failures[str(row["sample_id"])] = str(exc)
                print(
                    f"[{virus_id}] Phase 2 skipped: {sample_label} ({exc})",
                    file=sys.stderr,
                    flush=True,
                )

        if len(phase2_rows) < target_successes:
            print(
                f"[{virus_id}] Warning: only {len(phase2_rows)} of {target_successes} requested deep samples produced empirical TSS calls",
                file=sys.stderr,
                flush=True,
            )

    tss_values = [int(row["phase2_tss_coord"]) for row in phase2_rows]
    if not tss_values:
        failure_summary = "; ".join(
            f"{sample_id}: {error}" for sample_id, error in sorted(phase2_failures.items())
        )
        raise RuntimeError(
            f"No empirical TSS coordinates were detected for {virus_id}. Phase 2 failures: {failure_summary}"
        )
    consensus_tss = int(statistics.median(tss_values))

    combined_rows: list[dict[str, object]] = []
    selected_ids = {str(row["sample_id"]) for row in phase2_rows}
    phase2_by_sample = {str(row["sample_id"]): row for row in phase2_rows}
    for row in phase1_rows:
        merged = {
            "virus_id": virus_id,
            "dataset": row["dataset"],
            "sample_id": row["sample_id"],
            "phase1_pairs": row["phase1_pairs"],
            "phase1_total_reads": row["phase1_total_reads"],
            "phase1_mapped_primary_reads": row["phase1_mapped_primary_reads"],
            "phase1_viral_fraction": f"{float(row['phase1_viral_fraction']):.8f}",
            "selected_top_n": "yes" if str(row["sample_id"]) in selected_ids else "no",
            "phase2_status": "not_selected",
            "phase2_pairs": "",
            "phase2_tss_coord": "",
            "phase2_error": "",
        }
        if str(row["sample_id"]) in phase2_by_sample:
            merged["phase2_status"] = "selected"
            merged["phase2_pairs"] = phase2_by_sample[str(row["sample_id"])]["phase2_pairs"]
            merged["phase2_tss_coord"] = phase2_by_sample[str(row["sample_id"])]["phase2_tss_coord"]
        elif str(row["sample_id"]) in phase2_failures:
            merged["phase2_status"] = "failed"
            merged["phase2_error"] = phase2_failures[str(row["sample_id"])]
        combined_rows.append(merged)

    consensus_row = {
        "virus_id": virus_id,
        "source_fasta": ref_fasta.name,
        "source_cds_gtf": source_cds_gtf.name,
        "candidate_datasets": ",".join(spec["candidate_datasets"]),
        "selected_samples": ",".join(str(row["sample_id"]) for row in phase2_rows),
        "seqid": seqid,
        "genome_length": genome_length,
        "structural_start": structural_start,
        "search_start": search_start,
        "search_end": search_end,
        "phase1_pairs": args.phase1_pairs,
        "phase2_pairs": args.phase2_pairs,
        "top_n": args.top_n,
        "jump_window": args.jump_window,
        "jump_multiplier": args.jump_multiplier,
        "consensus_tss": consensus_tss,
    }

    curated_dir = output_root / virus_id
    curated_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(ref_fasta, curated_dir / "virus.fa")
    (curated_dir / "virus.gtf").write_text(
        build_curated_gtf(
            seqid=seqid,
            genome_length=genome_length,
            consensus_tss=consensus_tss,
            gene_id=args.gene_id or str(spec["gene_id"]),
            gene_name=str(spec["gene_name"]),
            tx_49s=str(spec["tx_49s"]),
            tx_26s=str(spec["tx_26s"]),
        ),
        encoding="utf-8",
    )

    return combined_rows, consensus_row


def main() -> int:
    args = parse_args()
    args.mapper_bin = require_command(args.mapper_bin)
    args.samtools_bin = require_command(args.samtools_bin)

    viruses = ("VEEV", "EEEV") if args.virus == "all" else (args.virus,)
    results_dir = Path(args.results_dir)
    all_rows: list[dict[str, object]] = []
    consensus_rows: list[dict[str, object]] = []

    for virus_id in viruses:
        rows, consensus = curate_one_virus(args, virus_id)
        all_rows.extend(rows)
        consensus_rows.append(consensus)

    write_results_tsv(results_dir / "empirical_tss_results.tsv", all_rows)
    write_consensus_tsv(results_dir / "consensus_tss.tsv", consensus_rows)

    for row in consensus_rows:
        print(
            f'{row["virus_id"]}: consensus 26S TSS {row["consensus_tss"]} '
            f'from {row["selected_samples"]}'
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
