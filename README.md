# V-EEEV Nat Hist nf-core/rnaseq ISAAC-NG Handoff

This scaffold runs `nf-core/rnaseq` version `3.23.0` for two separately sorted datasets:

- `mouse`
- `rat`

This pipeline expects the mixed local `V-EEEV Nat Hist` FASTQs to be sorted into species-specific input folders before upload to ISAAC-NG. The scripts do not infer species from filenames.

## Folder Contract

Populate the scaffold like this on ISAAC-NG before launching:

```text
veeev_nat_hist_nfcore/
├── bin/
├── inputs/
│   ├── mouse/
│   └── rat/
├── metadata/
├── settings.env
├── references/
│   ├── mouse/
│   │   ├── host/
│   │   ├── virus/
│   │   └── build/
│   └── rat/
│       ├── host/
│       ├── virus/
│       └── build/
├── results/
├── nextflow.config
└── submit_rnaseq.sh
```

Expected reference inputs for each species:

- `references/<species>/host/`: exactly one host FASTA and exactly one host GTF
- `references/<species>/virus/`: exactly one virus FASTA and exactly one virus GTF or GFF/GFF3

Generated files:

- `metadata/<species>_samplesheet.csv`
- `references/<species>/build/combined.fa`
- `references/<species>/build/combined.gtf`

Central settings:

- `settings.env`: account, module names, profile, scratch paths, and viral gene ID
- `submit_rnaseq.sh` `#SBATCH` lines: manager job partition, memory, cores, and walltime

## Upload Workflow

1. Sort FASTQs locally into `inputs/mouse/` and `inputs/rat/`.
2. Copy this entire folder to ISAAC-NG. Scripts may live in home or project storage, and large FASTQs may be staged wherever the run should read them from.
3. Drop the finalized host and virus references into the matching species reference folders.
4. Launch from a compute node through Slurm, not from the login node.

Example copy targets on ISAAC-NG:

- scripts and config: `/nfs/home/<netid>/veeev_nat_hist_nfcore`
- runtime output: `$SCRATCHDIR/veeev_nat_hist_nfcore`

## Reference Staging

The scaffold expects one host FASTA, one host GTF, one virus FASTA, and one virus annotation file per species. Host references are pinned here for reproducibility to the current Ensembl release used when this scaffold was prepared:

- mouse host: Ensembl release `115`, assembly `GRCm39`
- rat host: Ensembl release `115`, assembly `GRCr8`

The virus reference is still a placeholder and should be filled in once the final viral genome and annotation are selected.

Create the reference directories:

```bash
mkdir -p references/mouse/host references/mouse/virus
mkdir -p references/rat/host references/rat/virus
```

Pinned host download helper:

```bash
chmod +x bin/download_host_references.sh
bash bin/download_host_references.sh all
```

Virus download pattern with placeholder URLs:

```bash
curl -L "<mouse-virus-fasta-url>" -o references/mouse/virus/virus.fa
curl -L "<mouse-virus-gff-url>"   -o references/mouse/virus/virus.gff3

curl -L "<rat-virus-fasta-url>" -o references/rat/virus/virus.fa
curl -L "<rat-virus-gff-url>"   -o references/rat/virus/virus.gff3
```

The file basenames do not matter. The helper scripts only require that each folder contains exactly one matching file of the expected type.

ISAAC-NG transfer note:

- Download external references locally, not from the ISAAC login node.
- Upload references and FASTQs through the ISAAC data transfer nodes `dtn1.isaac.utk.edu` or `dtn2.isaac.utk.edu`.
- `dtn1` or `dtn2` should be used for `scp`, `sftp`, or `rsync`; `Globus` is preferred for large transfers.

Example upload to ISAAC using `rsync` through `dtn1`:

```bash
rsync -av references/mouse/host/ <netid>@dtn1.isaac.utk.edu:/nfs/home/<netid>/veeev-nat-hist-nfcore-isaac/references/mouse/host/
rsync -av references/rat/host/ <netid>@dtn1.isaac.utk.edu:/nfs/home/<netid>/veeev-nat-hist-nfcore-isaac/references/rat/host/
```

Reference preflight:

```bash
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh rat
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
export SHARED_VIRUS_GENE_ID="VEEV_SHARED_GENE"
```

`ISAAC_ACCOUNT` is required for Nextflow child jobs. `SBATCH_ACCOUNT` is recommended so the manager job also lands under the correct account when `sbatch submit_rnaseq.sh <species>` is used.

Notes:

- `SCRATCHDIR` must exist in the Slurm job environment. The scaffold no longer falls back to `/tmp`.
- `NFCORE_PROFILE=singularity` is the default unless the ISAAC environment requires something different.
- `DEFAULT_STRANDEDNESS` can be `auto`, `forward`, `reverse`, or `unstranded`.
- Leave `CONTAINER_MODULE` empty if `singularity` is already on `PATH` on ISAAC.
- The launcher checks `nextflow -version` and stops unless it is at least `25.04.3`, which is required by `nf-core/rnaseq` `3.23.0`.

## Run

Submit one species at a time:

```bash
sbatch submit_rnaseq.sh mouse
sbatch submit_rnaseq.sh rat
```

Optional lightweight preflight before submission:

```bash
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh mouse
PREFLIGHT_ONLY=1 bash submit_rnaseq.sh rat
```

That preflight checks:

- the species FASTQ folder exists and has `*_R1_001.fastq.gz` files
- the host and virus reference folders exist
- each reference folder has exactly one expected FASTA and annotation file

The launcher will:

1. build the combined host+virus reference for the selected species
2. auto-generate the nf-core samplesheet
3. run `nf-core/rnaseq` with `--aligner star_salmon`

## Quick Checks

Before launch:

```bash
ls -lh references/mouse/host references/mouse/virus
ls -lh references/rat/host references/rat/virus
```

After reference build:

```bash
ls -lh references/mouse/build/combined.fa references/mouse/build/combined.gtf
grep -c 'gene_id "VEEV_SHARED_GENE"' references/mouse/build/combined.gtf
```

Queue monitoring:

```bash
squeue -u <netid>
```

Dry-run the launcher logic without starting Nextflow:

```bash
DRY_RUN=1 SKIP_MODULE_LOAD=1 bash submit_rnaseq.sh mouse
```
