#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: make_samplesheet.sh <fastq_dir> <output_csv>

Scans a flat directory of paired-end FASTQs named:
  <sample>_R1_001.fastq.gz
  <sample>_R2_001.fastq.gz

Writes an nf-core/rnaseq samplesheet with columns:
  sample,fastq_1,fastq_2,strandedness

The script aborts on:
  - empty input folders
  - missing mates
  - duplicate sample stems
EOF
}

if [[ $# -ne 2 ]]; then
    usage >&2
    exit 1
fi

fastq_dir=$1
output_csv=$2
strandedness=${DEFAULT_STRANDEDNESS:-auto}

case "$strandedness" in
    auto|forward|reverse|unstranded) ;;
    *)
        echo "DEFAULT_STRANDEDNESS must be one of: auto, forward, reverse, unstranded" >&2
        exit 1
        ;;
esac

if [[ ! -d "$fastq_dir" ]]; then
    echo "FASTQ directory not found: $fastq_dir" >&2
    exit 1
fi

mkdir -p "$(dirname "$output_csv")"
tmp_output=$(mktemp "${output_csv}.XXXXXX")

mapfile -t r1_files < <(find "$fastq_dir" -maxdepth 1 \( -type f -o -type l \) -name '*_R1_001.fastq.gz' | sort)

if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "No R1 FASTQ files found in: $fastq_dir" >&2
    exit 1
fi

declare -A seen=()

{
    echo "sample,fastq_1,fastq_2,strandedness"

    for r1 in "${r1_files[@]}"; do
        base=$(basename "$r1")
        sample=${base%_R1_001.fastq.gz}

        if [[ -n ${seen["$sample"]+x} ]]; then
            echo "Duplicate sample stem detected: $sample" >&2
            exit 1
        fi
        seen["$sample"]=1

        r2="$fastq_dir/${sample}_R2_001.fastq.gz"
        if [[ ! -e "$r2" ]]; then
            echo "Missing mate for sample $sample: expected $r2" >&2
            exit 1
        fi

        printf '%s,%s,%s,%s\n' "$sample" "$r1" "$r2" "$strandedness"
    done
} >"$tmp_output"

expected_r2=${#r1_files[@]}
actual_r2=$(find "$fastq_dir" -maxdepth 1 \( -type f -o -type l \) -name '*_R2_001.fastq.gz' | wc -l | tr -d ' ')
if [[ "$actual_r2" -ne "$expected_r2" ]]; then
    rm -f "$tmp_output"
    echo "R2 file count ($actual_r2) does not match R1 file count ($expected_r2) in $fastq_dir" >&2
    exit 1
fi

mv "$tmp_output" "$output_csv"
echo "Wrote samplesheet: $output_csv"
