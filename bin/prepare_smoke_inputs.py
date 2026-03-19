#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create tiny paired FASTQ subsets for smoke-testing the real nf-core/rnaseq scaffold."
        )
    )
    parser.add_argument("dataset", choices=["mouse_veev", "mouse_eeev", "rat_veev"])
    parser.add_argument("source_dir", help="Flat source FASTQ directory for the selected dataset")
    parser.add_argument("output_dir", help="Destination flat FASTQ directory for smoke-test subsets")
    parser.add_argument(
        "--config",
        default=None,
        help="TSV manifest with columns: dataset, sample_id, read_pairs (default: metadata/smoke_samples.tsv)",
    )
    parser.add_argument(
        "--read-pairs",
        type=int,
        default=None,
        help="Override the default read-pair limit from the smoke manifest",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate smoke FASTQs even if the expected outputs already exist",
    )
    return parser.parse_args()


def load_manifest(path: Path) -> dict[str, tuple[str, int]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    required = {"dataset", "sample_id", "read_pairs"}
    if not rows or not required.issubset(reader.fieldnames or []):
        raise ValueError(f"Smoke manifest {path} must contain columns: dataset, sample_id, read_pairs")

    manifest: dict[str, tuple[str, int]] = {}
    for row in rows:
        dataset = row["dataset"].strip()
        if dataset in manifest:
            raise ValueError(f"Duplicate dataset entry in smoke manifest: {dataset}")
        manifest[dataset] = (row["sample_id"].strip(), int(row["read_pairs"]))
    return manifest


def resolve_default_manifest() -> Path:
    root = Path(__file__).resolve().parent.parent
    return root / "metadata" / "smoke_samples.tsv"


def expected_fastqs(sample_id: str, source_dir: Path) -> tuple[Path, Path]:
    return (
        source_dir / f"{sample_id}_R1_001.fastq.gz",
        source_dir / f"{sample_id}_R2_001.fastq.gz",
    )


def write_subset(r1_source: Path, r2_source: Path, r1_out: Path, r2_out: Path, read_pairs: int) -> int:
    tmp_r1 = r1_out.with_suffix(r1_out.suffix + ".tmp")
    tmp_r2 = r2_out.with_suffix(r2_out.suffix + ".tmp")
    written = 0

    try:
        with gzip.open(r1_source, "rt", encoding="utf-8") as r1_handle, gzip.open(
            r2_source, "rt", encoding="utf-8"
        ) as r2_handle, gzip.open(tmp_r1, "wt", encoding="utf-8") as r1_out_handle, gzip.open(
            tmp_r2, "wt", encoding="utf-8"
        ) as r2_out_handle:
            while written < read_pairs:
                r1_record = [r1_handle.readline() for _ in range(4)]
                r2_record = [r2_handle.readline() for _ in range(4)]

                r1_empty = all(line == "" for line in r1_record)
                r2_empty = all(line == "" for line in r2_record)
                if r1_empty and r2_empty:
                    break
                if r1_empty != r2_empty or any(line == "" for line in r1_record + r2_record):
                    raise ValueError("FASTQ mates ended unevenly while creating smoke subset")

                r1_out_handle.writelines(r1_record)
                r2_out_handle.writelines(r2_record)
                written += 1

        if written == 0:
            raise ValueError("No read pairs were written to the smoke subset")

        shutil.move(tmp_r1, r1_out)
        shutil.move(tmp_r2, r2_out)
        return written
    finally:
        for tmp_path in (tmp_r1, tmp_r2):
            if tmp_path.exists():
                tmp_path.unlink()


def main() -> int:
    args = parse_args()
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    config_path = Path(args.config) if args.config else resolve_default_manifest()

    if not source_dir.is_dir():
        raise FileNotFoundError(f"Smoke input source directory not found: {source_dir}")
    if not config_path.is_file():
        raise FileNotFoundError(f"Smoke manifest not found: {config_path}")

    manifest = load_manifest(config_path)
    if args.dataset not in manifest:
        raise ValueError(f"No smoke manifest entry found for dataset {args.dataset}")

    sample_id, default_read_pairs = manifest[args.dataset]
    read_pairs = args.read_pairs if args.read_pairs is not None else default_read_pairs
    if read_pairs <= 0:
        raise ValueError("read_pairs must be a positive integer")

    r1_source, r2_source = expected_fastqs(sample_id, source_dir)
    if not r1_source.is_file() or not r2_source.is_file():
        raise FileNotFoundError(
            f"Expected source FASTQs for smoke sample {sample_id} under {source_dir}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    r1_out = output_dir / r1_source.name
    r2_out = output_dir / r2_source.name

    if not args.force and r1_out.is_file() and r2_out.is_file():
        print(
            f"Smoke FASTQs already exist for {args.dataset} sample {sample_id}: "
            f"{r1_out} {r2_out}"
        )
        return 0

    written = write_subset(r1_source, r2_source, r1_out, r2_out, read_pairs)
    print(
        f"Created smoke FASTQs for {args.dataset} sample {sample_id} "
        f"with {written} read pairs in {output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
