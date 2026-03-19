# V-EEEV Nat Hist nf-core/rnaseq ISAAC-NG Handoff

This scaffold runs `nf-core/rnaseq` version `3.23.0` for the three dataset groups present in the `V-EEEV Nat Hist` delivery:

- `mouse_veev`
- `mouse_eeev`
- `rat_veev`

The mixed delivery can be classified automatically from the FASTQ stems. The staging command handles the source directory name `V-EEEV Nat Hist` safely as long as the path is quoted in the shell.

The viral references are data-driven and already tracked in this repo:

- raw viral source files live under `viral_references/raw/`
- empirical TSS ranking results live under `viral_references/curation/`
- curated pipeline-facing viral references live under `references/VEEV/` and `references/EEEV/`

## Folder Contract

Populate the scaffold like this on ISAAC-NG before launching:

```text
veeev_nat_hist_nfcore/
├── bin/
├── inputs/
│   ├── mouse_veev/
│   ├── mouse_eeev/
│   ├── rat_veev/
│   └── smoke/
│       ├── mouse_veev/
│       ├── mouse_eeev/
│       └── rat_veev/
├── metadata/
│   └── smoke/
├── settings.env
├── references/
│   ├── mouse/
│   ├── rat/
│   ├── EEEV/
│   ├── VEEV/
│   └── build/
├── viral_references/
│   ├── raw/
│   └── curation/
├── results/
├── nextflow.config
└── submit_rnaseq.sh
```

Expected reference inputs:

- `references/mouse/`: exactly one mouse host FASTA and exactly one mouse host GTF
- `references/rat/`: exactly one rat host FASTA and exactly one rat host GTF
- `references/EEEV/`: curated empirical `virus.fa` and `virus.gtf`
- `references/VEEV/`: curated empirical `virus.fa` and `virus.gtf`

Generated files:

- `metadata/nat_hist_manifest.csv`
- `metadata/catalog/nat_hist_samples.tsv`
- `metadata/catalog/<dataset>.tsv`
- `metadata/<dataset>_samplesheet.csv`
- `metadata/smoke/<dataset>_samplesheet.csv`
- `references/build/<dataset>/combined.fa`
- `references/build/<dataset>/combined.gtf`
- `viral_references/curation/empirical_tss_results.tsv`
- `viral_references/curation/consensus_tss.tsv`

Central settings:

- `settings.env`: account, module names, profile, scratch paths, and viral gene ID
- `submit_rnaseq.sh` `#SBATCH` lines: manager job partition, memory, cores, and walltime

## Upload Workflow

1. Copy this entire folder to ISAAC-NG.
2. Stage the mixed delivery into `inputs/mouse_veev/`, `inputs/mouse_eeev/`, and `inputs/rat_veev/` with `bin/stage_nat_hist_inputs.py`.
3. Drop the finalized host and virus references into the matching reference folders.
4. Launch from a compute node through Slurm, not from the login node.

Example copy targets on ISAAC-NG:

- scripts and config: `/nfs/home/<netid>/veeev_nat_hist_nfcore`
- runtime output: `$SCRATCHDIR/veeev_nat_hist_nfcore`

## Reference Staging

The scaffold expects separate host and virus reference folders. Host references are pinned here for reproducibility to the current Ensembl release used when this scaffold was prepared:

- mouse host: Ensembl release `115`, assembly `GRCm39`
- rat host: Ensembl release `115`, assembly `GRCr8`

The host references are downloaded into `references/mouse/` and `references/rat/`. The viral references are already tracked in the repo as raw source files plus curated empirical outputs.

Create the reference directories:

```bash
mkdir -p references/mouse references/rat references/EEEV references/VEEV
```

Pinned host download helper:

```bash
chmod +x bin/download_host_references.sh
bash bin/download_host_references.sh all
```

Virus download pattern with placeholder URLs:

```bash
ls -lh viral_references/raw/EEEV
ls -lh viral_references/raw/VEEV
```

The committed raw viral source files are:

- `viral_references/raw/VEEV/VEEV_TrD.fa`
- `viral_references/raw/VEEV/VEEV_TrD (CDS).gtf`
- `viral_references/raw/EEEV/EEEV_FL93.fa`
- `viral_references/raw/EEEV/EEEV_FL93 (CDS).gtf`

The committed curated viral outputs are:

- `references/VEEV/virus.fa`
- `references/VEEV/virus.gtf`
- `references/EEEV/virus.fa`
- `references/EEEV/virus.gtf`

These curated viral GTFs contain exactly two transcript models per virus:

- full-length genomic RNA `49S`
- empirical subgenomic RNA `26S`

Both transcripts share one `gene_id` and keep distinct `transcript_id` values.

Current empirical consensus TSS values:

- `VEEV_26S` starts at `7567`
- `EEEV_26S` starts at `7550`

ISAAC-NG transfer note:

- Download external references locally, not from the ISAAC login node.
- Upload references and FASTQs through the ISAAC data transfer nodes `dtn1.isaac.utk.edu` or `dtn2.isaac.utk.edu`.
- `dtn1` or `dtn2` should be used for `scp`, `sftp`, or `rsync`; `Globus` is preferred for large transfers.

Example upload to ISAAC using `rsync` through `dtn1`:

```bash
rsync -av references/mouse/ <netid>@dtn1.isaac.utk.edu:/nfs/home/<netid>/veeev-nat-hist-nfcore-isaac/references/mouse/
rsync -av references/rat/ <netid>@dtn1.isaac.utk.edu:/nfs/home/<netid>/veeev-nat-hist-nfcore-isaac/references/rat/
rsync -av references/EEEV/ <netid>@dtn1.isaac.utk.edu:/nfs/home/<netid>/veeev-nat-hist-nfcore-isaac/references/EEEV/
rsync -av references/VEEV/ <netid>@dtn1.isaac.utk.edu:/nfs/home/<netid>/veeev-nat-hist-nfcore-isaac/references/VEEV/
rsync -av viral_references/ <netid>@dtn1.isaac.utk.edu:/nfs/home/<netid>/veeev-nat-hist-nfcore-isaac/viral_references/
```

## Viral Curation

The viral references are not taken directly from the source GTF/GFF coordinates. They are derived from the RNA-seq data with a two-tier viral-only pilot mapping workflow implemented in `bin/curate_viral_references.py`.

Phase 1:

- map the first `100k` read pairs from every candidate sample to the virus only
- rank samples by `viral mapped reads / total sampled reads`

Phase 2:

- deep-map the first `1M` read pairs from the top `5` ranked samples
- calculate `samtools depth` across the structural-junction search window
- detect the first strong coverage jump as the empirical `26S` TSS
- take the median of the successful top-sample calls as the final consensus

Virus pools:

- `VEEV`: `mouse_veev` + `rat_veev`
- `EEEV`: `mouse_eeev`

Required local tools for curation:

- `minimap2`
- `samtools`
- `python3`

Example rerun command from the repo root:

```bash
python3 bin/curate_viral_references.py all \
  --mapper-bin /home/alex_ubuntu/miniconda3/envs/vpipe/bin/minimap2 \
  --samtools-bin /home/alex_ubuntu/miniconda3/envs/vpipe/bin/samtools \
  --jobs 4
```

Tracked outputs written by that command:

- `viral_references/curation/empirical_tss_results.tsv`
- `viral_references/curation/consensus_tss.tsv`
- `references/VEEV/virus.fa`
- `references/VEEV/virus.gtf`
- `references/EEEV/virus.fa`
- `references/EEEV/virus.gtf`

Input staging from the mixed delivery:

```bash
python3 bin/stage_nat_hist_inputs.py "/path/to/V-EEEV Nat Hist"
```

Filename rules used by the staging command:

- `100-199` and `300-399` -> `mouse_veev`
- `200-299`, `400-499`, and `4xxxx` -> `mouse_eeev`
- `B*` and `L*` -> `rat_veev`

That command will:

- write `metadata/nat_hist_manifest.csv`
- classify the mixed FASTQs into `mouse_veev`, `mouse_eeev`, and `rat_veev`
- stage symlinks under `inputs/` without duplicating the original FASTQs

Optional staging flags:

```bash
python3 bin/stage_nat_hist_inputs.py "/path/to/V-EEEV Nat Hist" --manifest-only
python3 bin/stage_nat_hist_inputs.py "/path/to/V-EEEV Nat Hist" --method copy --clean
```

## Metadata Catalog

The scaffold now keeps two separate metadata products for the natural history inputs:

- `metadata/nat_hist_manifest.csv`: local staging manifest with machine-specific absolute paths; this file is intentionally ignored by Git
- `metadata/catalog/*.tsv`: tracked, repo-safe sample catalog with stable annotations and repo-relative staged FASTQ paths

Export the tracked catalog from the local manifest with:

```bash
python3 bin/export_nat_hist_catalog.py
```

Optional explicit paths:

```bash
python3 bin/export_nat_hist_catalog.py \
  --manifest metadata/nat_hist_manifest.csv \
  --output-dir metadata/catalog
```

Tracked catalog outputs:

- `metadata/catalog/nat_hist_samples.tsv`
- `metadata/catalog/mouse_veev.tsv`
- `metadata/catalog/mouse_eeev.tsv`
- `metadata/catalog/rat_veev.tsv`

Column notes:

- `study_code`: mouse study block code when it is explicitly known
- `route`: inoculation route for mapped mouse study blocks
- `day`: intentionally blank until an explicit per-sample day mapping exists
- `infection_status`: `unknown` for mouse studies until the virus/mock map is resolved; blank for the separate rat study

Current metadata defaults:

- three-digit mouse samples in studies `045`, `046`, `047`, and `048` keep blank `day` and `infection_status=unknown`
- five-digit `4xxxx` samples are treated as a separate `EEEV BDGR251` study, not as part of study `048`
- `B*` and `L*` samples are tracked as a separate `B/L Chibuke Rat study`

Reference preflight:

```bash
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_veev
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_eeev
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh rat_veev
```

## Smoke Testing

Use smoke mode to run one tiny representative sample per dataset through the full `nf-core/rnaseq` workflow before launching the real jobs. The smoke inputs are generated on demand from the staged real FASTQs and written under `inputs/smoke/`.

Tracked smoke configuration:

- `metadata/smoke_samples.tsv`: representative sample ID and default read-pair count per dataset
- `metadata/smoke_container_urls.txt`: known container URLs to pre-cache on the login node

Representative smoke samples currently configured:

- `mouse_veev` -> `101`
- `mouse_eeev` -> `201`
- `rat_veev` -> `B10`

Generate or refresh one smoke subset and verify preflight:

```bash
SMOKE_TEST=1 PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_veev
SMOKE_TEST=1 PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_eeev
SMOKE_TEST=1 PREFLIGHT_ONLY=1 bash submit_rnaseq.sh rat_veev
```

Override the default subset size or force regeneration:

```bash
SMOKE_TEST=1 SMOKE_READ_PAIRS=50000 SMOKE_REGENERATE=1 PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_veev
```

Smoke runs use separate scratch-backed output roots:

- `results_smoke/<dataset>`
- `work_smoke/<dataset>`

Run all three smoke tests sequentially:

```bash
job1=$(sbatch --parsable --account="$ISAAC_ACCOUNT" --export=ALL,SMOKE_TEST=1 submit_rnaseq.sh mouse_veev)
job2=$(sbatch --parsable --account="$ISAAC_ACCOUNT" --dependency=afterany:$job1 --export=ALL,SMOKE_TEST=1 submit_rnaseq.sh mouse_eeev)
sbatch --account="$ISAAC_ACCOUNT" --dependency=afterany:$job2 --export=ALL,SMOKE_TEST=1 submit_rnaseq.sh rat_veev
```

A passing smoke test should exercise the full workflow, including:

- `SALMON_QUANT`
- `STAR_ALIGN`
- `SUBREAD_FEATURECOUNTS`
- `UCSC_BEDCLIP` / `BEDGRAPH...BIGWIG`

## Container Pre-cache

On clusters where compute nodes cannot reach the internet, pre-pull known container images from the login node before smoke or production launches:

```bash
bash bin/precache_nfcore_containers.sh
```

To augment the known URL list with images seen in prior manager logs:

```bash
bash bin/precache_nfcore_containers.sh \
  --log-file nfcore_rnaseq.5069375.out \
  --log-file nfcore_rnaseq.5080364.out
```

## Environment

`settings.env` is the main place to change runtime settings.

The manager job resources near the top of `submit_rnaseq.sh` control the `#SBATCH` settings.

Shell overrides are still supported, but the default workflow keeps normal values in `settings.env`.

Example settings:

```bash
export ISAAC_ACCOUNT="ACF-UTKXXXX"
export SBATCH_ACCOUNT="$ISAAC_ACCOUNT"
export NEXTFLOW_MODULE="nextflow"
export CONTAINER_MODULE=""
export NFCORE_PROFILE="singularity"
export DEFAULT_STRANDEDNESS="auto"
export RESULTS_BASE="$SCRATCHDIR/veeev_nat_hist_nfcore/results"
export WORK_ROOT="$SCRATCHDIR/veeev_nat_hist_nfcore/work"
export CONTAINER_CACHE="$SCRATCHDIR/veeev_nat_hist_nfcore/containers"
export NXF_HOME="$SCRATCHDIR/veeev_nat_hist_nfcore/.nextflow"
```

`ISAAC_ACCOUNT` is required for Nextflow child jobs. `SBATCH_ACCOUNT` is recommended so the manager job also lands under the correct account when `sbatch submit_rnaseq.sh <dataset>` is used.

Notes:

- `SCRATCHDIR` must exist in the Slurm job environment. The scaffold no longer falls back to `/tmp`.
- `NFCORE_PROFILE=singularity` is the default unless the ISAAC environment requires something different.
- `DEFAULT_STRANDEDNESS` can be `auto`, `forward`, `reverse`, or `unstranded`.
- Leave `CONTAINER_MODULE` empty if `singularity` is already on `PATH` on ISAAC.
- The launcher checks `nextflow -version` and stops unless it is at least `25.04.3`, which is required by `nf-core/rnaseq` `3.23.0`.

## Run

Submit one dataset at a time:

```bash
sbatch submit_rnaseq.sh mouse_veev
sbatch submit_rnaseq.sh mouse_eeev
sbatch submit_rnaseq.sh rat_veev
```

Optional lightweight preflight before submission:

```bash
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_veev
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse_eeev
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh rat_veev
```

That preflight checks:

- the dataset FASTQ folder exists and has `*_R1_001.fastq.gz` files
- the host and virus reference folders exist
- each reference folder has exactly one expected FASTA and annotation file

The launcher will:

1. build the combined host+virus reference for the selected dataset
2. auto-generate the nf-core samplesheet
3. run `nf-core/rnaseq` with `--aligner star_salmon`

## Quick Checks

Before launch:

```bash
ls -lh references/mouse references/rat references/EEEV references/VEEV
```

After reference build:

```bash
ls -lh references/build/mouse_veev/combined.fa references/build/mouse_veev/combined.gtf
grep -c 'gene_id "VEEV"' references/build/mouse_veev/combined.gtf
grep -c 'gene_id "EEEV"' references/build/mouse_eeev/combined.gtf
grep -c 'transcript_id "VEEV_26S"' references/build/mouse_veev/combined.gtf
grep -c 'transcript_id "EEEV_26S"' references/build/mouse_eeev/combined.gtf
```

Queue monitoring:

```bash
squeue -u <netid>
```

Dry-run the launcher logic without starting Nextflow:

```bash
DRY_RUN=1 SKIP_MODULE_LOAD=1 bash submit_rnaseq.sh mouse_veev
```
