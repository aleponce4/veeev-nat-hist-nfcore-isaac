#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: download_host_references.sh [mouse|rat|all]

Downloads pinned Ensembl host references into:
  references/mouse/
  references/rat/

Default: all

Pinned versions:
  mouse -> Ensembl release 115, GRCm39
  rat   -> Ensembl release 115, GRCr8
EOF
}

target=${1:-all}
case "$target" in
    mouse|rat|all) ;;
    -h|--help)
        usage
        exit 0
        ;;
    *)
        usage >&2
        exit 1
        ;;
esac

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd "$script_dir/.." && pwd)

download_file() {
    local url=$1
    local dest=$2
    mkdir -p "$(dirname "$dest")"
    printf 'Downloading %s -> %s\n' "$url" "$dest"
    curl -fL "$url" -o "$dest"
}

download_file_with_fallback() {
    local dest=$1
    shift
    local url
    for url in "$@"; do
        printf 'Trying %s -> %s\n' "$url" "$dest"
        if curl -fL "$url" -o "$dest"; then
            return 0
        fi
    done
    echo "All download URLs failed for $dest" >&2
    exit 1
}

download_mouse() {
    download_file \
        "https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" \
        "$root_dir/references/mouse/mouse.fa.gz"
    download_file \
        "https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz" \
        "$root_dir/references/mouse/mouse.gtf.gz"
}

download_rat() {
    download_file_with_fallback \
        "$root_dir/references/rat/rat.fa.gz" \
        "https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.primary_assembly.fa.gz" \
        "https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz"
    download_file \
        "https://ftp.ensembl.org/pub/release-115/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.115.gtf.gz" \
        "$root_dir/references/rat/rat.gtf.gz"
}

case "$target" in
    mouse) download_mouse ;;
    rat) download_rat ;;
    all)
        download_mouse
        download_rat
        ;;
esac

printf 'Host reference download complete for: %s\n' "$target"
