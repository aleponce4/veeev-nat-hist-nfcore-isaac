#!/usr/bin/env bash
#SBATCH -J nfcore_rnaseq
#SBATCH -p campus
#SBATCH -q campus
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH -o %x.%j.out

set -euo pipefail

# Main places to edit settings:
# 1. The #SBATCH lines above control the manager job submitted by sbatch.
# 2. settings.env controls the reusable runtime settings used below.
# 3. CONTAINER_MODULE can be left empty if ISAAC already provides singularity on PATH.
#
# Common HPC terminal commands:
#   cd /path/to/veeev_nat_hist_nfcore
#   nano settings.env
#   module avail nextflow
#   python3 bin/stage_nat_hist_inputs.py "/path/to/V-EEEV Nat Hist"
#   PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_veev
#   export ISAAC_ACCOUNT="ACF-UTKXXXX"
#   export SBATCH_ACCOUNT="$ISAAC_ACCOUNT"
#   DRY_RUN=1 SKIP_MODULE_LOAD=1 bash submit_rnaseq.sh mouse_veev
#   sbatch submit_rnaseq.sh mouse_veev
#   sbatch submit_rnaseq.sh mouse_eeev
#   sbatch submit_rnaseq.sh rat_veev
#   squeue -u <netid>

usage() {
    cat <<'EOF'
Usage:
  sbatch submit_rnaseq.sh <mouse_veev|mouse_eeev|rat_veev>

Primary settings:
  edit settings.env

Required final value:
  ISAAC_ACCOUNT must be changed from the template account

Quick HPC commands:
  cd /path/to/veeev_nat_hist_nfcore
  nano settings.env
  python3 bin/stage_nat_hist_inputs.py "/path/to/V-EEEV Nat Hist"
  PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_veev
  export ISAAC_ACCOUNT="ACF-UTKXXXX"
  export SBATCH_ACCOUNT="$ISAAC_ACCOUNT"
  sbatch submit_rnaseq.sh mouse_veev
EOF
}

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 1
fi

# Step 1: identify which dataset to run.
dataset=$1
case "$dataset" in
    mouse_veev|mouse_eeev|rat_veev) ;;
    *)
        echo "Dataset must be one of: mouse_veev, mouse_eeev, rat_veev. Got: $dataset" >&2
        exit 1
        ;;
esac

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
settings_file="$script_dir/settings.env"
preflight_only=${PREFLIGHT_ONLY:-0}

# Step 2: load the central runtime settings, if present.
if [[ -f "$settings_file" ]]; then
    # shellcheck disable=SC1091
    source "$settings_file"
fi

# Step 3: define the dataset-specific inputs, references, and output targets.
samplesheet="$script_dir/metadata/${dataset}_samplesheet.csv"
fastq_dir="$script_dir/inputs/$dataset"
reference_builder="$script_dir/bin/build_combined_reference.sh"
samplesheet_builder="$script_dir/bin/make_samplesheet.sh"
case "$dataset" in
    mouse_veev)
        host_ref="mouse"
        virus_ref="VEEV"
        ;;
    mouse_eeev)
        host_ref="mouse"
        virus_ref="EEEV"
        ;;
    rat_veev)
        host_ref="rat"
        virus_ref="VEEV"
        ;;
esac
combined_fasta="$script_dir/references/build/$dataset/combined.fa"
combined_gtf="$script_dir/references/build/$dataset/combined.gtf"
host_dir="$script_dir/references/$host_ref"
virus_dir="$script_dir/references/$virus_ref"

nextflow_module=${NEXTFLOW_MODULE:-}
container_module=${CONTAINER_MODULE:-}
nfcore_profile=${NFCORE_PROFILE:-}
results_base=${RESULTS_BASE:-}
work_root=${WORK_ROOT:-}
container_cache=${CONTAINER_CACHE:-}

work_dir="$work_root/$dataset"
outdir="$results_base/$dataset"

find_single_matching_file() {
    local search_dir=$1
    shift
    local patterns=("$@")
    local matches=()
    local pattern

    for pattern in "${patterns[@]}"; do
        while IFS= read -r file; do
            matches+=("$file")
        done < <(find "$search_dir" -maxdepth 1 -type f -name "$pattern" | sort)
    done

    if [[ ${#matches[@]} -eq 0 ]]; then
        echo "No matching files found in $search_dir for patterns: ${patterns[*]}" >&2
        exit 1
    fi
    if [[ ${#matches[@]} -ne 1 ]]; then
        echo "Expected exactly one matching file in $search_dir, found ${#matches[@]}:" >&2
        printf '  %s\n' "${matches[@]}" >&2
        exit 1
    fi
}

check_fastq_inputs() {
    local first_r1

    if [[ ! -d "$fastq_dir" ]]; then
        echo "FASTQ input directory not found: $fastq_dir" >&2
        exit 1
    fi
    first_r1=$(find "$fastq_dir" -maxdepth 1 \( -type f -o -type l \) -name '*_R1_001.fastq.gz' -print -quit)
    if [[ -z "$first_r1" ]]; then
        echo "No R1 FASTQ files found in: $fastq_dir" >&2
        exit 1
    fi
}

check_reference_inputs() {
    if [[ ! -d "$host_dir" || ! -d "$virus_dir" ]]; then
        echo "Expected reference directories are missing for $dataset" >&2
        echo "Expected: $host_dir and $virus_dir" >&2
        exit 1
    fi

    find_single_matching_file "$host_dir" '*.fa' '*.fa.gz' '*.fasta' '*.fasta.gz' '*.fna' '*.fna.gz'
    find_single_matching_file "$host_dir" '*.gtf' '*.gtf.gz'
    find_single_matching_file "$virus_dir" '*.fa' '*.fa.gz' '*.fasta' '*.fasta.gz' '*.fna' '*.fna.gz'
    find_single_matching_file "$virus_dir" '*.gtf' '*.gtf.gz' '*.gff' '*.gff.gz' '*.gff3' '*.gff3.gz'
}

# Step 4: run lightweight preflight checks that can be used before sbatch.
check_fastq_inputs
check_reference_inputs

if [[ "$preflight_only" == "1" ]]; then
    printf 'Preflight checks passed for %s\n' "$dataset"
    printf 'FASTQs: %s\n' "$fastq_dir"
    printf 'Host reference dir:  %s (%s)\n' "$host_dir" "$host_ref"
    printf 'Virus reference dir: %s (%s)\n' "$virus_dir" "$virus_ref"
    exit 0
fi

# Step 5: validate the cluster environment needed for a real run.
if [[ -z "${ISAAC_ACCOUNT:-}" ]]; then
    echo "ISAAC_ACCOUNT is required" >&2
    exit 1
fi

if [[ "${ISAAC_ACCOUNT}" == "ACF-UTKXXXX" ]]; then
    echo "ISAAC_ACCOUNT is still set to the template value in settings.env" >&2
    exit 1
fi

if [[ -z "${SCRATCHDIR:-}" ]]; then
    echo "SCRATCHDIR is required on ISAAC-NG" >&2
    exit 1
fi

if [[ -z "$results_base" || -z "$work_root" || -z "$container_cache" ]]; then
    echo "Scratch-backed runtime paths are empty; confirm SCRATCHDIR is exported before launch" >&2
    exit 1
fi

export NXF_HOME=${NXF_HOME:-$SCRATCHDIR/veeev_nat_hist_nfcore/.nextflow}
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-$container_cache}
export NXF_APPTAINER_CACHEDIR=${NXF_APPTAINER_CACHEDIR:-$container_cache}

# Step 6: create the directories Nextflow and the helper scripts will write to.
mkdir -p "$results_base" "$work_dir" "$container_cache" "$script_dir/metadata"

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

nextflow_version_ge() {
    local current=$1
    local minimum=$2
    local current_major current_minor current_patch
    local minimum_major minimum_minor minimum_patch

    IFS=. read -r current_major current_minor current_patch <<<"$current"
    IFS=. read -r minimum_major minimum_minor minimum_patch <<<"$minimum"

    current_patch=${current_patch:-0}
    minimum_patch=${minimum_patch:-0}

    if (( current_major > minimum_major )); then
        return 0
    fi
    if (( current_major < minimum_major )); then
        return 1
    fi
    if (( current_minor > minimum_minor )); then
        return 0
    fi
    if (( current_minor < minimum_minor )); then
        return 1
    fi
    (( current_patch >= minimum_patch ))
}

check_nextflow_version() {
    local required_minimum="25.04.3"
    local detected_version

    if ! command -v nextflow >/dev/null 2>&1; then
        echo "nextflow is not on PATH after module loading" >&2
        exit 1
    fi

    detected_version=$(nextflow -version 2>/dev/null | awk '/version/ {print $NF; exit}')
    detected_version=${detected_version#v}

    if [[ -z "$detected_version" ]]; then
        echo "Could not determine Nextflow version from 'nextflow -version'" >&2
        exit 1
    fi

    if ! nextflow_version_ge "$detected_version" "$required_minimum"; then
        echo "Nextflow $detected_version is too old for nf-core/rnaseq 3.23.0; need >= $required_minimum" >&2
        exit 1
    fi
}

# Step 7: load Nextflow and the container runtime unless the environment already has them.
if [[ "${SKIP_MODULE_LOAD:-0}" != "1" ]]; then
    if ensure_module_command; then
        module purge
        module load "$nextflow_module"
        if [[ -n "$container_module" ]]; then
            module load "$container_module"
        fi
    else
        echo "Environment module command not found; set SKIP_MODULE_LOAD=1 only if Nextflow is already on PATH" >&2
        exit 1
    fi
fi

# Step 7b: fail early if the loaded Nextflow is too old for this nf-core release.
check_nextflow_version

# Step 8: build the combined reference and generate the nf-core samplesheet.
"$reference_builder" "$dataset"
"$samplesheet_builder" "$fastq_dir" "$samplesheet"

# Step 9: assemble the nf-core/rnaseq command with dataset-specific paths.
cmd=(
    nextflow run nf-core/rnaseq
    -r 3.23.0
    -profile "$nfcore_profile"
    -c "$script_dir/nextflow.config"
    -work-dir "$work_dir"
    --input "$samplesheet"
    --fasta "$combined_fasta"
    --gtf "$combined_gtf"
    --aligner star_salmon
    --outdir "$outdir"
    -resume
)

# Optional debugging mode to print the exact Nextflow command and stop.
if [[ "${DRY_RUN:-0}" == "1" ]]; then
    printf 'DRY_RUN nextflow command:\n'
    printf '%q ' "${cmd[@]}"
    printf '\n'
    exit 0
fi

# Step 10: launch the actual pipeline on the allocated compute node.
printf 'Launching nf-core/rnaseq for %s\n' "$dataset"
printf 'Results: %s\n' "$outdir"
printf 'Work dir: %s\n' "$work_dir"
printf 'Settings file: %s\n' "$settings_file"

"${cmd[@]}"
