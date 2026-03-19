#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  precache_nfcore_containers.sh [--cache-dir DIR] [--url-file FILE] [--log-file FILE ...] [--print-only]

Pre-pulls known nf-core/rnaseq container images into the shared cache from a login node.
Container URLs are read from:
  - a tracked URL list file
  - optional nfcore_rnaseq.*.out logs that contain "Pulling Singularity image ..."
EOF
}

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd "$script_dir/.." && pwd)
cache_dir=${CONTAINER_CACHE:-}
if [[ -z "$cache_dir" && -n "${SCRATCHDIR:-}" ]]; then
    cache_dir="${SCRATCHDIR}/veeev_nat_hist_nfcore/containers"
fi
url_file="$root_dir/metadata/smoke_container_urls.txt"
print_only=0
log_files=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cache-dir)
            cache_dir=$2
            shift 2
            ;;
        --url-file)
            url_file=$2
            shift 2
            ;;
        --log-file)
            log_files+=("$2")
            shift 2
            ;;
        --print-only)
            print_only=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "${cache_dir:-}" ]]; then
    echo "Container cache directory is empty; set CONTAINER_CACHE, SCRATCHDIR, or --cache-dir" >&2
    exit 1
fi

declare -A urls=()

if [[ -f "$url_file" ]]; then
    while IFS= read -r line; do
        [[ -n "${line// }" ]] || continue
        [[ $line == \#* ]] && continue
        urls["$line"]=1
    done <"$url_file"
fi

for log_file in "${log_files[@]}"; do
    [[ -f "$log_file" ]] || continue
    while IFS= read -r line; do
        if [[ $line =~ Pulling\ Singularity\ image\ (https://[^[:space:]]+) ]]; then
            urls["${BASH_REMATCH[1]}"]=1
        fi
    done <"$log_file"
done

if [[ ${#urls[@]} -eq 0 ]]; then
    echo "No container URLs were discovered; provide --url-file and/or --log-file inputs" >&2
    exit 1
fi

runtime=
if [[ $print_only -ne 1 ]]; then
    if command -v singularity >/dev/null 2>&1; then
        runtime=singularity
    elif command -v apptainer >/dev/null 2>&1; then
        runtime=apptainer
    else
        echo "Neither singularity nor apptainer is available on PATH" >&2
        exit 1
    fi
fi

mkdir -p "$cache_dir"
cd "$cache_dir"

image_name_from_url() {
    local url=$1
    local stripped=${url#https://}
    stripped=${stripped#http://}
    stripped=${stripped//\//-}
    stripped=${stripped//:/-}
    printf '%s.img\n' "$stripped"
}

for url in "${!urls[@]}"; do
    image_name=$(image_name_from_url "$url")
    if [[ $print_only -eq 1 ]]; then
        printf '%s\t%s\n' "$image_name" "$url"
        continue
    fi
    if [[ -f "$image_name" ]]; then
        printf 'Already cached: %s\n' "$image_name"
        continue
    fi
    printf 'Pulling %s -> %s\n' "$url" "$image_name"
    "$runtime" pull --name "$image_name" "$url"
done
