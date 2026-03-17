#!/usr/bin/env python3
"""Classify and stage the mixed V-EEEV Nat Hist delivery into dataset inputs."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys
from collections import Counter
from pathlib import Path


FASTQ_RE = re.compile(r"(?P<sample>.+)_R(?P<read>[12])_001\.fastq\.gz$")


def classify_sample(sample_id: str) -> tuple[str, str, str, str]:
    if re.fullmatch(r"[BL]\d+", sample_id):
        return ("rat", "VEEV", "rat_veev", "rat_letter_code")

    if re.fullmatch(r"\d{5}", sample_id) and sample_id.startswith("4"):
        return ("mouse", "EEEV", "mouse_eeev", "five_digit_4xxxx")

    if re.fullmatch(r"\d{3}", sample_id):
        numeric_id = int(sample_id)
        if 100 <= numeric_id <= 199:
            return ("mouse", "VEEV", "mouse_veev", "three_digit_100_199")
        if 200 <= numeric_id <= 299:
            return ("mouse", "EEEV", "mouse_eeev", "three_digit_200_299")
        if 300 <= numeric_id <= 399:
            return ("mouse", "VEEV", "mouse_veev", "three_digit_300_399")
        if 400 <= numeric_id <= 499:
            return ("mouse", "EEEV", "mouse_eeev", "three_digit_400_499")

    raise ValueError(f"Unrecognized sample naming pattern: {sample_id}")


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description=(
            "Scan a mixed V-EEEV Nat Hist FASTQ delivery, write a manifest, and "
            "stage dataset-specific inputs for the nf-core scaffold."
        )
    )
    parser.add_argument(
        "source_dir",
        help='Path to the mixed delivery folder, e.g. "/path/to/V-EEEV Nat Hist"',
    )
    parser.add_argument(
        "--output-root",
        default=str(repo_root),
        help="Scaffold root that contains inputs/ and metadata/ (default: repository root)",
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="Manifest CSV path (default: <output-root>/metadata/nat_hist_manifest.csv)",
    )
    parser.add_argument(
        "--method",
        choices=("symlink", "hardlink", "copy"),
        default="symlink",
        help="How to stage FASTQs into inputs/ (default: symlink)",
    )
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="Write the manifest but do not stage FASTQs into inputs/",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove existing staged FASTQs in the destination inputs/ folders before restaging",
    )
    return parser.parse_args()


def scan_delivery(source_dir: Path) -> list[dict[str, str]]:
    samples: dict[str, dict[str, Path]] = {}

    for path in sorted(source_dir.iterdir()):
        if not path.is_file():
            continue
        match = FASTQ_RE.fullmatch(path.name)
        if not match:
            continue

        sample_id = match.group("sample")
        read = match.group("read")
        sample_reads = samples.setdefault(sample_id, {})
        if read in sample_reads:
            raise RuntimeError(f"Duplicate FASTQ detected for sample {sample_id} read {read}: {path}")
        sample_reads[read] = path.resolve()

    if not samples:
        raise RuntimeError(f"No FASTQ files matching *_R1_001.fastq.gz were found in {source_dir}")

    manifest_rows = []
    for sample_id in sorted(samples):
        reads = samples[sample_id]
        if set(reads) != {"1", "2"}:
            present = ",".join(f"R{read}" for read in sorted(reads))
            raise RuntimeError(f"Incomplete read pair for sample {sample_id}: found {present}")

        host_ref, virus_ref, dataset, rule = classify_sample(sample_id)
        manifest_rows.append(
            {
                "sample_id": sample_id,
                "host_ref": host_ref,
                "virus_ref": virus_ref,
                "dataset": dataset,
                "classification_rule": rule,
                "source_fastq_1": str(reads["1"]),
                "source_fastq_2": str(reads["2"]),
            }
        )

    return manifest_rows


def ensure_clean_destinations(planned_targets: dict[Path, Path], clean: bool) -> None:
    dataset_dirs = sorted({path.parent for path in planned_targets})

    for dataset_dir in dataset_dirs:
        dataset_dir.mkdir(parents=True, exist_ok=True)
        existing_fastqs = {
            path.name: path
            for path in dataset_dir.iterdir()
            if path.is_file() and path.name.endswith(".fastq.gz")
        }
        expected_fastqs = {
            path.name: path for path in planned_targets if path.parent == dataset_dir
        }

        unexpected = sorted(set(existing_fastqs) - set(expected_fastqs))
        conflicting = []
        for name in sorted(set(existing_fastqs) & set(expected_fastqs)):
            existing = existing_fastqs[name]
            target = expected_fastqs[name]
            source = planned_targets[target]
            if existing.is_symlink():
                current_target = Path(os.path.realpath(existing))
                if current_target == source:
                    continue
            else:
                try:
                    if os.path.samefile(existing, source):
                        continue
                except FileNotFoundError:
                    pass
            conflicting.append(name)

        if not clean and (unexpected or conflicting):
            details = []
            if unexpected:
                details.append(f"unexpected files: {', '.join(unexpected)}")
            if conflicting:
                details.append(f"conflicting files: {', '.join(conflicting)}")
            raise RuntimeError(
                f"Existing staged FASTQs found in {dataset_dir}. "
                f"Rerun with --clean to replace them ({'; '.join(details)})."
            )

        if clean:
            for name in sorted(set(existing_fastqs)):
                existing_fastqs[name].unlink()


def stage_one_file(source: Path, destination: Path, method: str) -> None:
    if destination.exists() or destination.is_symlink():
        return

    if method == "symlink":
        destination.symlink_to(source)
    elif method == "hardlink":
        os.link(source, destination)
    elif method == "copy":
        shutil.copy2(source, destination)
    else:
        raise ValueError(f"Unsupported staging method: {method}")


def write_manifest(manifest_path: Path, rows: list[dict[str, str]]) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "host_ref",
        "virus_ref",
        "dataset",
        "classification_rule",
        "source_fastq_1",
        "source_fastq_2",
        "staged_fastq_1",
        "staged_fastq_2",
    ]
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()

    source_dir = Path(args.source_dir).expanduser().resolve()
    output_root = Path(args.output_root).expanduser().resolve()
    manifest_path = (
        Path(args.manifest).expanduser().resolve()
        if args.manifest
        else output_root / "metadata" / "nat_hist_manifest.csv"
    )

    if not source_dir.is_dir():
        raise RuntimeError(f"Source directory not found: {source_dir}")

    manifest_rows = scan_delivery(source_dir)

    planned_targets: dict[Path, Path] = {}
    for row in manifest_rows:
        dataset_dir = output_root / "inputs" / row["dataset"]
        target_r1 = dataset_dir / Path(row["source_fastq_1"]).name
        target_r2 = dataset_dir / Path(row["source_fastq_2"]).name
        row["staged_fastq_1"] = str(target_r1)
        row["staged_fastq_2"] = str(target_r2)
        planned_targets[target_r1] = Path(row["source_fastq_1"])
        planned_targets[target_r2] = Path(row["source_fastq_2"])

    if not args.manifest_only:
        ensure_clean_destinations(planned_targets, clean=args.clean)
        for destination, source in sorted(planned_targets.items()):
            destination.parent.mkdir(parents=True, exist_ok=True)
            stage_one_file(source, destination, args.method)

    write_manifest(manifest_path, manifest_rows)

    counts = Counter(row["dataset"] for row in manifest_rows)
    print(f"Manifest: {manifest_path}")
    if args.manifest_only:
        print("Staging: skipped (--manifest-only)")
    else:
        print(f"Staging method: {args.method}")
    for dataset in sorted(counts):
        print(f"{dataset}: {counts[dataset]} samples")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
