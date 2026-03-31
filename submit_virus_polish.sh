#!/usr/bin/env bash
#SBATCH -J virus_polish
#SBATCH -p campus
#SBATCH -q campus
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH -o %x.%j.out

set -euo pipefail

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

usage() {
    cat <<'EOF'
Usage:
  sbatch submit_virus_polish.sh <VEEV|EEEV|all>

Optional environment overrides:
  POLISH_RESULTS_ROOT        existing nf-core results root
  POLISH_OUTPUT_ROOT         output root for polished FASTAs and summaries
  POLISH_TOP_N              top-ranked BAMs to pool per virus
  POLISH_JOBS               concurrent worker count inside the polish script
  POLISH_MIN_MAPPED_READS   minimum viral mapped reads to keep a BAM
  POLISH_MASK_DEPTH         mask consensus positions below this depth with N
  POLISH_MAX_DEPTH          bcftools mpileup max depth cap
  INSTALL_POLISHED_REFERENCE=1  overwrite references/<virus>/virus.fa after backup
EOF
}

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 1
fi

virus=$1
case "$virus" in
    VEEV|EEEV|all) ;;
    *)
        echo "Virus must be one of: VEEV, EEEV, all. Got: $virus" >&2
        exit 1
        ;;
esac

if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/settings.env" ]]; then
    script_dir=$(cd "$SLURM_SUBMIT_DIR" && pwd)
else
    script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fi

settings_file="$script_dir/settings.env"
if [[ -f "$settings_file" ]]; then
    # shellcheck disable=SC1091
    source "$settings_file"
fi

results_root=${POLISH_RESULTS_ROOT:-${RESULTS_BASE:-$script_dir/results}}
output_root=${POLISH_OUTPUT_ROOT:-$script_dir/viral_references/polish}
samtools_module=${SAMTOOLS_MODULE:-samtools}
bcftools_module=${BCFTOOLS_MODULE:-bcftools}

ensure_module_command() {
    if command -v module >/dev/null 2>&1; then
        return 0
    fi
    if [[ -f /etc/profile.d/modules.sh ]]; then
        # shellcheck disable=SC1091
        source /etc/profile.d/modules.sh
    elif [[ -f /etc/profile.d/lmod.sh ]]; then
        # shellcheck disable=SC1091
        source /etc/profile.d/lmod.sh
    fi
    command -v module >/dev/null 2>&1
}

if [[ "${SKIP_MODULE_LOAD:-0}" != "1" ]]; then
    if ensure_module_command; then
        module purge
        module load "$samtools_module"
        module load "$bcftools_module"
    else
        echo "Environment module command not found; set SKIP_MODULE_LOAD=1 only if samtools and bcftools are already on PATH" >&2
        exit 1
    fi
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "samtools is not on PATH after module loading" >&2
    exit 1
fi
if ! command -v bcftools >/dev/null 2>&1; then
    echo "bcftools is not on PATH after module loading" >&2
    exit 1
fi

cmd=(
    python3 "$script_dir/bin/polish_virus_reference_from_bams.py"
    "$virus"
    --results-root "$results_root"
    --reference-root "$script_dir/references"
    --output-root "$output_root"
    --top-n "${POLISH_TOP_N:-5}"
    --jobs "${POLISH_JOBS:-${SLURM_CPUS_PER_TASK:-16}}"
    --min-mapped-reads "${POLISH_MIN_MAPPED_READS:-0}"
    --mask-depth "${POLISH_MASK_DEPTH:-10}"
    --max-depth "${POLISH_MAX_DEPTH:-1000000}"
)

if [[ "${INSTALL_POLISHED_REFERENCE:-0}" == "1" ]]; then
    cmd+=(--install)
fi

if [[ "${DRY_RUN:-0}" == "1" ]]; then
    printf 'DRY_RUN polish command:\n'
    printf '%q ' "${cmd[@]}"
    printf '\n'
    exit 0
fi

printf 'Polishing viral reference(s): %s\n' "$virus"
printf 'Results root: %s\n' "$results_root"
printf 'Output root: %s\n' "$output_root"
printf 'Settings file: %s\n' "$settings_file"

"${cmd[@]}"
