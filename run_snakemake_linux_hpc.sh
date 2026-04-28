#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-/Users/evanelko/Documents/pepseq_snakemake}"
SNAKEFILE="${SNAKEFILE:-$WORKDIR/Snakefile}"
CONFIGFILE="${CONFIGFILE:-$WORKDIR/config.yaml}"
PROFILE="${SNAKEMAKE_PROFILE:-}"
JOBS="${SNAKEMAKE_JOBS:-${SLURM_CPUS_PER_TASK:-32}}"
LATENCY_WAIT="${SNAKEMAKE_LATENCY_WAIT:-120}"
RESTART_TIMES="${SNAKEMAKE_RESTART_TIMES:-1}"

if ! command -v conda >/dev/null 2>&1; then
  echo "Error: conda not found in PATH." >&2
  exit 1
fi

if ! command -v snakemake >/dev/null 2>&1; then
  echo "Error: snakemake not found in PATH." >&2
  exit 1
fi

# Force conda environment resolution/creation for Linux HPC nodes.
export CONDA_SUBDIR=linux-64

echo "Using CONDA_SUBDIR=${CONDA_SUBDIR}"
echo "Running Snakemake in: ${WORKDIR}"
echo "Using configfile: ${CONFIGFILE}"
echo "Jobs: ${JOBS}"
echo "Latency wait: ${LATENCY_WAIT}s"

cd "$WORKDIR"

if [[ -n "$PROFILE" ]]; then
  echo "Using profile: ${PROFILE}"
fi

if [[ -n "$PROFILE" && "$PROFILE" == *"slurm"* ]]; then
  if ! snakemake --help 2>&1 | grep -Eq -- "--executor .*\\{[^}]*slurm"; then
    echo "Error: Slurm profile selected, but this Snakemake install has no Slurm executor plugin." >&2
    echo "Install one of:" >&2
    echo "  pip install snakemake-executor-plugin-slurm" >&2
    echo "  mamba install -c conda-forge -c bioconda snakemake-executor-plugin-slurm" >&2
    exit 1
  fi
fi

CMD=(
  snakemake
  --use-conda
  --rerun-incomplete
  --keep-going
  --printshellcmds
  --jobs "$JOBS"
  --latency-wait "$LATENCY_WAIT"
  --restart-times "$RESTART_TIMES"
  --replace-workflow-config
  -s "$SNAKEFILE"
  --configfile "$CONFIGFILE"
)

if [[ -n "$PROFILE" ]]; then
  CMD+=(--profile "$PROFILE")
fi

"${CMD[@]}" "$@"
