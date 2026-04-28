# PepSeq Snakemake Pipeline

## Overview
This repository runs a QC + binning pipeline for PepSeq count matrices, with optional AutoPepsirf and PSEA stages.

Main workflow file:
- `Snakefile`

Main config file:
- `config.yaml`

Launcher scripts:
- `run_snakemake_osx64.sh` (local/mac launcher)
- `run_snakemake_linux_hpc.sh` (Linux/HPC launcher)

Profiles:
- `profiles/local/config.yaml`
- `profiles/slurm/config.yaml`

## Pipeline Stages
1. Shorten sample column names
2. Read-count threshold filtering
3. Zero-count distribution filtering (with optional histogram)
4. Pairwise Pearson correlation filtering
5. QC failure summary merge
6. AutomateBins control selection + bin generation
7. Optional AutoPepsirf run
8. Optional PSEA preparation + run

## Python Scripts (Pipeline)
- `python_scripts/shorten_column_names.py`: shorten sample names and write mapping table
- `python_scripts/QC_RC_thresh.py`: read-count threshold QC + pass/fail outputs
- `python_scripts/downsample2max_zero_filter.py`: downsample, zero-count QC, and optional all-sample/downsampled outputs
- `python_scripts/pairwise_pearson_from_pairs.py`: pairwise Pearson QC from pairs table
- `python_scripts/merge_qc_failures.py`: merge per-step QC failure logs
- `python_scripts/filter_pairs_by_columns.py`: subset pairs to samples present in a matrix
- `python_scripts/remove_control_samples.py`: remove selected control columns and restore full names
- `python_scripts/automateBins.py`: choose control Sblks and generate bins

## Requirements
At minimum, the environment running Snakemake and the Python scripts should include:
- `snakemake`
- `python`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `snakemake-executor-plugin-slurm` (for SLURM profile)
- `matrixtools` (from LadnerLab Modules)

Also required:
- `pepsirf` CLI available to rules/scripts that call it
- `conda` in `PATH`

## Quick Start
### Local dry-run
```bash
/Users/evanelko/Documents/pepseq_snakemake/run_snakemake_osx64.sh --dry-run
```

### Full DAG dry-run (repo-local validation)
Use this when your default `config.yaml` points to external files that are not present locally.
```bash
snakemake --dry-run --snakefile Snakefile --configfile config.dryrun.yaml
```

### Local run
```bash
/Users/evanelko/Documents/pepseq_snakemake/run_snakemake_osx64.sh
```

### HPC dry-run (SLURM profile)
```bash
SNAKEMAKE_PROFILE=/path/to/snakemake/pepsirf_snakemake_pipeline/profiles/slurm \
/path/to/snakemake/pepsirf_snakemake_pipeline/run_snakemake_linux_hpc.sh --dry-run
```

### HPC run (SLURM profile)
```bash
SNAKEMAKE_PROFILE=/path/to/snakemake/pepsirf_snakemake_pipeline/profiles/slurm \
/path/to/snakemake/pepsirf_snakemake_pipeline/run_snakemake_linux_hpc.sh --rerun-triggers mtime
```

## Running With `tmux` on HPC (recommended)
```bash
ssh <cluster>
tmux new -s snakemake
cd /path/to/snakemake/pepsirf_snakemake_pipeline
SNAKEMAKE_PROFILE=/path/to/snakemake/pepsirf_snakemake_pipeline/profiles/slurm \
./run_snakemake_linux_hpc.sh --rerun-triggers mtime
```

Detach:
```bash
Ctrl-b d
```

Reattach:
```bash
tmux attach -t snakemake
```

## Configuration Reference (`config.yaml`)
### Required inputs
- `raw_counts_tsv`: raw count matrix TSV input
- `pairs_tsv`: pair definitions TSV
- `sblk_id`: prefix used to identify superblock controls (example: `Sblk`)

### QC thresholds and toggles
- `rc_cutoff`: minimum readcount threshold for RC filter
- `run_correlation_qc`: enable/disable pairwise correlation QC
- `pearson_cutoff`: Pearson cutoff for correlation QC
- `run_zscore_correlation_qc`: optional Z-score correlation QC (requires AutoPepsirf)
- `zscore_pearson_cutoff`: Pearson cutoff for Z-score QC
- `force_rerun_with_token`: if `true`, injects run-token dependency to force recomputation
- `downsample_max`: max reads used in downsampling
- `downsample_min`: min reads allowed before sample exclusion
- `zero_sd_mult`: SD multiplier for zero-count cutoff
- `zero_filter_seed`: random seed for deterministic downsampling
- `zero_hist_output`: histogram output path; supports `"{output_dir}/..."` placeholder

Zero-count QC stage outputs:
- `02_zero_filtered.tsv`: matrix after zero-count filtering
- `02_zero_failed_samples.tsv`: samples that failed zero-count filtering
- `02_zero_all_counts.tsv`: zero-count values for all samples considered by this step (pass + fail)
- `02_zero_downsampled_counts.tsv`: full downsampled counts matrix used for zero-count QC (before pass/fail filtering)

### Correlation scatterplot options
- `correlation_scatterplots`: enable scatterplot generation
- `correlation_scatterplot_dir`: output dir for plots
- `correlation_plot_format`: `png`, `pdf`, or `svg`
- `correlation_plot_dpi`: figure DPI
- `correlation_max_pairs`: max plotted pairs (`null` for all)

### Script path options
- `python_scripts_dir`: default directory for pipeline Python scripts
- Optional overrides (commented by default):
  - `rc_script`
  - `correlation_script`
  - `zero_filter_script`
  - `automate_bins_script`
  - `qc_failure_merge_script`
  - `shorten_columns_script`
  - `remove_controls_script`
  - `filter_pairs_script`

### Per-rule resources
- `rule_resources`: mapping of rule name to:
  - `threads`
  - `mem_mb`
  - `runtime` (minutes)

### Output controls
- `output_dir`: root directory for `qc`, `logs`, `run`, `autopepsirf`, `psea`
- `automatebins_outprefix`: output prefix for generated bin files

### AutomateBins options
- `automatebins_top_fraction`
- `automatebins_min_readcount_sum`
- `automatebins_no_scatterplots`
- `automatebins_max_pairs`
- `automatebins_plot_format`
- `automatebins_plot_dpi`

### AutoPepsirf options (used when `run_autopepsirf: true`)
- `run_autopepsirf`
- `autopepsirf_conda_path`
- `autopepsirf_negative_id`
- `autopepsirf_exact_z_thresh`
- `autopepsirf_exact_cs_thresh`
- `autopepsirf_raw_constraint`
- `autopepsirf_hdi`
- `autopepsirf_exact_zenrich_thresh`
- `autopepsirf_tsv_base_str`
- `autopepsirf_pepsirf_tsv_dir`

### PSEA options (used when `run_psea: true`)
- `run_psea` (requires `run_autopepsirf: true`)
- `psea_conda_path`
- `psea_epitope_metadata`
- `psea_peptide_sets_gmt`
- `psea_species_tsv`
- `psea_threshold`
- `psea_pval_thresh`
- `psea_nes_thresh`
- `psea_min_size`
- `psea_max_size`
- `psea_permutation_num`
- `psea_spline_type`
- `psea_degree`
- `psea_vis_outputs_dir`
- `psea_iterative_analysis`
- `psea_seed`

## Profile Reference
### Local profile (`profiles/local/config.yaml`)
- `cores: 4`
- `jobs: 4`
- `use-conda: true`
- restart/latency/keep-going defaults for local execution

### SLURM profile (`profiles/slurm/config.yaml`)
- `executor: slurm`
- `jobs: 100`
- `use-conda: true`
- `default-resources: mem_mb=4000, runtime=60`
- Optional cluster-specific settings (recommended):
  - `slurm-account`
  - `slurm-partition`
  - `slurm-qos`

## Launcher Script Environment Variables
`run_snakemake_linux_hpc.sh` supports:
- `WORKDIR`
- `SNAKEFILE`
- `CONFIGFILE`
- `SNAKEMAKE_PROFILE`
- `SNAKEMAKE_JOBS`
- `SNAKEMAKE_LATENCY_WAIT`
- `SNAKEMAKE_RESTART_TIMES`

Example:
```bash
WORKDIR=/path/to/snakemake/pepsirf_snakemake_pipeline \
SNAKEFILE=/path/to/snakemake/pepsirf_snakemake_pipeline/Snakefile \
CONFIGFILE=/path/to/snakemake/pepsirf_snakemake_pipeline/config_hpc.yaml \
SNAKEMAKE_PROFILE=/path/to/snakemake/pepsirf_snakemake_pipeline/profiles/slurm \
SNAKEMAKE_JOBS=100 \
/path/to/snakemake/pepsirf_snakemake_pipeline/run_snakemake_linux_hpc.sh --dry-run
```

## Common Commands
### Unlock after interrupted run
```bash
/path/to/snakemake/pepsirf_snakemake_pipeline/run_snakemake_linux_hpc.sh --unlock
```

### Re-run from incomplete jobs
`run_snakemake_linux_hpc.sh` already uses `--rerun-incomplete` by default.

### Force one setting at runtime
```bash
/path/to/snakemake/pepsirf_snakemake_pipeline/run_snakemake_linux_hpc.sh \
--config raw_counts_tsv=/path/to/data/table.tsv
```

## Troubleshooting
### `ModuleNotFoundError: seaborn`
Install in the Snakemake runtime environment:
```bash
conda install -n snakemake -c conda-forge seaborn -y
```

### `invalid choice: 'slurm'`
Install SLURM executor plugin:
```bash
conda install -n snakemake -c conda-forge -c bioconda snakemake-executor-plugin-slurm -y
```

### `Directory cannot be locked`
Run unlock once, then rerun:
```bash
./run_snakemake_linux_hpc.sh --unlock
```

### Missing input files
Confirm `raw_counts_tsv` and `pairs_tsv` in config are valid absolute paths on the host where Snakemake runs.

## Pre-Push Checklist
Before pushing this repository to GitHub:
1. Run a dry run:
```bash
/Users/evanelko/Documents/pepseq_snakemake/run_snakemake_osx64.sh --dry-run
```
2. Optionally run the pipeline for a full integration check:
```bash
/Users/evanelko/Documents/pepseq_snakemake/run_snakemake_osx64.sh
```
3. Confirm generated outputs and logs are ignored (`results/`, `.snakemake/`, daily `*.log`) via `.gitignore`.
4. Verify `config.yaml` does not contain machine-specific secrets/paths you do not want in version control.
