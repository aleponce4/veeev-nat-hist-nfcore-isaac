# V-EEEV Nat Hist RNA-seq Pipeline

Small wrapper repo for running `nf-core/rnaseq` (`3.23.0`) on the V-EEEV natural history datasets on HPC.

It does four things:

1. stages a mixed FASTQ delivery into dataset-specific folders
2. builds combined host + virus references
3. generates nf-core samplesheets
4. launches `nf-core/rnaseq` for one of these datasets:
   - `mouse_veev`
   - `mouse_eeev`
   - `rat_veev`

## Expected Input

FASTQs:

- Input FASTQs are paired-end and named like `SAMPLE_R1_001.fastq.gz` and `SAMPLE_R2_001.fastq.gz`
- If starting from the mixed `V-EEEV Nat Hist` delivery, use `bin/stage_nat_hist_inputs.py`
- That script classifies samples into:
  - `mouse_veev`
  - `mouse_eeev`
  - `rat_veev`

References:

- `references/mouse/`: one mouse FASTA and one mouse GTF
- `references/rat/`: one rat FASTA and one rat GTF
- `references/VEEV/`: `virus.fa` and `virus.gtf`
- `references/EEEV/`: `virus.fa` and `virus.gtf`

The repo already includes curated viral references in [references/VEEV/virus.fa](/mnt/c/Users/alepo/Documents/veeev-nat-hist-nfcore-isaac/references/VEEV/virus.fa), [references/VEEV/virus.gtf](/mnt/c/Users/alepo/Documents/veeev-nat-hist-nfcore-isaac/references/VEEV/virus.gtf), [references/EEEV/virus.fa](/mnt/c/Users/alepo/Documents/veeev-nat-hist-nfcore-isaac/references/EEEV/virus.fa), and [references/EEEV/virus.gtf](/mnt/c/Users/alepo/Documents/veeev-nat-hist-nfcore-isaac/references/EEEV/virus.gtf).

## Basic Workflow

Stage the mixed delivery:

```bash
python3 bin/stage_nat_hist_inputs.py "/path/to/V-EEEV Nat Hist"
```

Download host references if needed:

```bash
bash bin/download_host_references.sh all
```

Edit cluster/runtime settings:

```bash
nano settings.env
```

Important setting:

- set `ISAAC_ACCOUNT` to the real account instead of `ACF-UTKXXXX`

Preflight check a dataset:

```bash
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_veev
```

Run a dataset:

```bash
sbatch submit_rnaseq.sh mouse_veev
sbatch submit_rnaseq.sh mouse_eeev
sbatch submit_rnaseq.sh rat_veev
```

## What `submit_rnaseq.sh` Does

For the selected dataset, the launcher:

1. checks inputs and references
2. builds `references/build/<dataset>/combined.fa`
3. builds `references/build/<dataset>/combined.gtf`
4. writes `metadata/<dataset>_samplesheet.csv`
5. runs `nextflow run nf-core/rnaseq -r 3.23.0`

The nf-core run uses:

- `--aligner star_salmon`
- `-profile "$NFCORE_PROFILE"`
- `-c nextflow.config`

## Useful Commands

Make a samplesheet manually:

```bash
bash bin/make_samplesheet.sh inputs/mouse_veev metadata/mouse_veev_samplesheet.csv
```

Build a combined reference manually:

```bash
bash bin/build_combined_reference.sh mouse_veev
```

Print the exact Nextflow command without running:

```bash
DRY_RUN=1 SKIP_MODULE_LOAD=1 bash submit_rnaseq.sh mouse_veev
```

Run a smoke test:

```bash
SMOKE_TEST=1 PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_veev
SMOKE_TEST=1 sbatch submit_rnaseq.sh mouse_veev
```

Pre-cache nf-core containers:

```bash
bash bin/precache_nfcore_containers.sh
```

## Output Locations

Main outputs go to scratch-backed paths from `settings.env`:

- `RESULTS_BASE/<dataset>`
- `WORK_ROOT/<dataset>`
- `CONTAINER_CACHE`

Smoke test outputs go to:

- `RESULTS_BASE_SMOKE/<dataset>`
- `WORK_ROOT_SMOKE/<dataset>`

## Requirements

- Slurm
- Nextflow `>= 25.04.3`
- Singularity or Apptainer
- Python 3

This repo is meant as a practical run wrapper, not a general-purpose distributed pipeline package.
