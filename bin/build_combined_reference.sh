#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: build_combined_reference.sh <mouse|rat>

Builds:
  references/<species>/build/combined.fa
  references/<species>/build/combined.gtf

Expected inputs:
  references/<species>/host/   -> exactly one FASTA and exactly one GTF
  references/<species>/virus/  -> exactly one FASTA and exactly one GTF/GFF/GFF3
EOF
}

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 1
fi

species=$1
case "$species" in
    mouse|rat) ;;
    *)
        echo "Species must be 'mouse' or 'rat', got: $species" >&2
        exit 1
        ;;
esac

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd "$script_dir/.." && pwd)
species_dir="$root_dir/references/$species"
host_dir="$species_dir/host"
virus_dir="$species_dir/virus"
build_dir="$species_dir/build"
shared_gene_id=${SHARED_VIRUS_GENE_ID:-VEEV_SHARED_GENE}

if [[ ! -d "$host_dir" || ! -d "$virus_dir" ]]; then
    echo "Expected reference directories are missing under $species_dir" >&2
    exit 1
fi

find_single_file() {
    local search_dir=$1
    shift
    local patterns=("$@")
    local files=()
    local pattern
    for pattern in "${patterns[@]}"; do
        while IFS= read -r file; do
            files+=("$file")
        done < <(find "$search_dir" -maxdepth 1 -type f -name "$pattern" | sort)
    done

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "No matching files found in $search_dir for patterns: ${patterns[*]}" >&2
        exit 1
    fi
    if [[ ${#files[@]} -ne 1 ]]; then
        echo "Expected exactly one matching file in $search_dir, found ${#files[@]}:" >&2
        printf '  %s\n' "${files[@]}" >&2
        exit 1
    fi
    printf '%s\n' "${files[0]}"
}

stream_file() {
    local path=$1
    if [[ $path == *.gz ]]; then
        gzip -cd "$path"
    else
        cat "$path"
    fi
}

require_command() {
    local cmd=$1
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Required command not found on PATH: $cmd" >&2
        exit 1
    fi
}

check_duplicate_fasta_ids() {
    local duplicates
    duplicates=$(
        {
            stream_file "$host_fasta"
            stream_file "$virus_fasta"
        } | awk '
            /^>/ {
                id = substr($1, 2)
                if (++seen[id] == 2) {
                    print id
                }
            }
        '
    )

    if [[ -n "$duplicates" ]]; then
        echo "Duplicate FASTA sequence IDs found across host and virus references:" >&2
        printf '  %s\n' $duplicates >&2
        exit 1
    fi
}

require_command python3

host_fasta=$(find_single_file "$host_dir" '*.fa' '*.fa.gz' '*.fasta' '*.fasta.gz' '*.fna' '*.fna.gz')
host_gtf=$(find_single_file "$host_dir" '*.gtf' '*.gtf.gz')
virus_fasta=$(find_single_file "$virus_dir" '*.fa' '*.fa.gz' '*.fasta' '*.fasta.gz' '*.fna' '*.fna.gz')
virus_annotation=$(find_single_file "$virus_dir" '*.gtf' '*.gtf.gz' '*.gff' '*.gff.gz' '*.gff3' '*.gff3.gz')

check_duplicate_fasta_ids

mkdir -p "$build_dir"

combined_fasta="$build_dir/combined.fa"
combined_gtf="$build_dir/combined.gtf"
normalized_virus_gtf="$build_dir/virus.normalized.gtf"

{
    stream_file "$host_fasta"
    printf '\n'
    stream_file "$virus_fasta"
    printf '\n'
} >"$combined_fasta"

python3 "$script_dir/normalize_virus_gtf.py" \
    "$virus_annotation" \
    "$normalized_virus_gtf" \
    --shared-gene-id "$shared_gene_id"

if [[ ! -s "$normalized_virus_gtf" ]]; then
    echo "Normalized viral GTF is missing or empty: $normalized_virus_gtf" >&2
    exit 1
fi

{
    stream_file "$host_gtf"
    printf '\n'
    cat "$normalized_virus_gtf"
} >"$combined_gtf"

if [[ ! -s "$combined_fasta" || ! -s "$combined_gtf" ]]; then
    echo "Combined reference build failed for $species" >&2
    exit 1
fi

echo "Built combined FASTA: $combined_fasta"
echo "Built combined GTF:   $combined_gtf"
