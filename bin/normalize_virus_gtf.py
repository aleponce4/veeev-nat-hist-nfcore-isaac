#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path


TRANSCRIPT_FEATURES = {"transcript", "mrna"}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Normalize a viral GTF/GFF so all viral transcripts share one gene_id "
            "while keeping distinct transcript_id values."
        )
    )
    parser.add_argument("input_path", help="Input virus annotation in GTF, GFF, or GFF3 format")
    parser.add_argument("output_path", help="Output GTF path")
    parser.add_argument(
        "--shared-gene-id",
        default="VEEV_SHARED_GENE",
        help='gene_id to assign to all viral features (default: "VEEV_SHARED_GENE")',
    )
    return parser.parse_args()


def parse_gtf_attributes(raw: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for chunk in raw.strip().strip(";").split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if " " in chunk:
            key, value = chunk.split(" ", 1)
            attributes[key] = value.strip().strip('"')
        elif "=" in chunk:
            key, value = chunk.split("=", 1)
            attributes[key] = value.strip().strip('"')
    return attributes


def parse_gff_attributes(raw: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for chunk in raw.strip().split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" in chunk:
            key, value = chunk.split("=", 1)
            attributes[key] = value
    return attributes


def format_gtf_attributes(attributes: dict[str, str]) -> str:
    ordered = []
    for key, value in attributes.items():
        escaped = value.replace("\\", "\\\\").replace('"', '\\"')
        ordered.append(f'{key} "{escaped}"')
    return "; ".join(ordered) + ";"


def choose_value(attributes: dict[str, str], *keys: str) -> str | None:
    for key in keys:
        value = attributes.get(key)
        if value:
            return value
    return None


def sanitize_identifier(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", value).strip("_") or "viral_feature"


def detect_format(lines: list[str]) -> str:
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        attrs = fields[8]
        if "gene_id" in attrs and '"' in attrs:
            return "gtf"
        if "=" in attrs:
            return "gff"
    raise ValueError("Could not determine whether the virus annotation is GTF or GFF/GFF3")


def normalize_gtf(lines: list[str], shared_gene_id: str) -> list[str]:
    normalized: list[str] = []
    for index, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        if line.startswith("#"):
            normalized.append(line if line.endswith("\n") else line + "\n")
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"Expected 9 tab-delimited columns in GTF line {index}")

        feature = fields[2].lower()
        attrs = parse_gtf_attributes(fields[8])
        gene_name = choose_value(attrs, "gene_name", "gene", "Name", "locus_tag", "gene_id")
        transcript_id = choose_value(
            attrs,
            "transcript_id",
            "transcript_name",
            "Name",
            "ID",
            "Parent",
            "locus_tag",
            "gene_id",
        )

        out_attrs: dict[str, str] = {"gene_id": shared_gene_id}
        if gene_name:
            out_attrs["gene_name"] = gene_name

        if feature in TRANSCRIPT_FEATURES:
            fields[2] = "transcript"

        if feature != "gene":
            if not transcript_id:
                transcript_id = f"{fields[2]}_{index}"
            transcript_id = sanitize_identifier(transcript_id)
            out_attrs["transcript_id"] = transcript_id
            out_attrs["transcript_name"] = choose_value(attrs, "transcript_name", "Name") or transcript_id

        fields[8] = format_gtf_attributes(out_attrs)
        normalized.append("\t".join(fields) + "\n")

    return normalized


def normalize_gff(lines: list[str], shared_gene_id: str) -> list[str]:
    records = []
    id_to_tx: dict[str, str] = {}
    id_to_gene_name: dict[str, str] = {}
    tx_to_gene_ref: dict[str, str] = {}
    tx_to_tx_name: dict[str, str] = {}

    for index, line in enumerate(lines, start=1):
        if not line.strip() or line.startswith("#"):
            records.append({"raw": line, "comment": True})
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"Expected 9 tab-delimited columns in GFF line {index}")

        attrs = parse_gff_attributes(fields[8])
        feature = fields[2].lower()
        record = {"fields": fields, "attrs": attrs, "feature": feature, "index": index, "comment": False}
        records.append(record)

        feature_id = choose_value(attrs, "ID")
        if feature == "gene" and feature_id:
            id_to_gene_name[feature_id] = choose_value(attrs, "Name", "gene", "gene_name") or feature_id

        if feature in TRANSCRIPT_FEATURES:
            transcript_id = sanitize_identifier(
                choose_value(attrs, "transcript_id", "ID", "Name", "transcript", "Parent") or f"transcript_{index}"
            )
            if feature_id:
                id_to_tx[feature_id] = transcript_id
            if parent := choose_value(attrs, "Parent"):
                tx_to_gene_ref[transcript_id] = parent.split(",")[0]
            tx_to_tx_name[transcript_id] = choose_value(attrs, "Name") or transcript_id

    normalized: list[str] = []
    for record in records:
        if record["comment"]:
            raw = record["raw"]
            normalized.append(raw if raw.endswith("\n") else raw + "\n")
            continue

        fields = list(record["fields"])
        attrs = record["attrs"]
        feature = record["feature"]
        index = record["index"]

        parent_values = [item for item in attrs.get("Parent", "").split(",") if item]
        transcript_id = None
        for parent in parent_values:
            if parent in id_to_tx:
                transcript_id = id_to_tx[parent]
                break

        if feature in TRANSCRIPT_FEATURES:
            transcript_id = sanitize_identifier(
                choose_value(attrs, "transcript_id", "ID", "Name", "transcript") or f"transcript_{index}"
            )
            fields[2] = "transcript"
        elif feature != "gene" and not transcript_id:
            transcript_id = sanitize_identifier(
                choose_value(attrs, "transcript_id", "Name", "ID", "Parent") or f"{fields[2]}_{index}"
            )

        gene_name = None
        for parent in parent_values:
            if parent in id_to_gene_name:
                gene_name = id_to_gene_name[parent]
                break
        if transcript_id and transcript_id in tx_to_gene_ref and not gene_name:
            gene_name = id_to_gene_name.get(tx_to_gene_ref[transcript_id])
        gene_name = gene_name or choose_value(attrs, "gene", "gene_name", "Name")

        out_attrs: dict[str, str] = {"gene_id": shared_gene_id}
        if gene_name:
            out_attrs["gene_name"] = gene_name
        if feature != "gene":
            out_attrs["transcript_id"] = transcript_id or f"{fields[2]}_{index}"
            out_attrs["transcript_name"] = (
                tx_to_tx_name.get(out_attrs["transcript_id"])
                or choose_value(attrs, "Name")
                or out_attrs["transcript_id"]
            )

        fields[8] = format_gtf_attributes(out_attrs)
        normalized.append("\t".join(fields) + "\n")

    return normalized


def main() -> int:
    args = parse_args()
    input_path = Path(args.input_path)
    output_path = Path(args.output_path)

    with open_text(input_path) as handle:
        lines = handle.readlines()

    file_format = detect_format(lines)
    if file_format == "gtf":
        normalized = normalize_gtf(lines, args.shared_gene_id)
    else:
        normalized = normalize_gff(lines, args.shared_gene_id)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.writelines(normalized)

    transcript_ids = set()
    with output_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            match = re.search(r'transcript_id "([^"]+)"', line)
            if match:
                transcript_ids.add(match.group(1))

    if not transcript_ids:
        raise ValueError("No transcript_id values were generated in the normalized viral GTF")

    print(
        f"Normalized viral annotation to {output_path} with shared gene_id "
        f'"{args.shared_gene_id}" and {len(transcript_ids)} transcript(s).',
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
