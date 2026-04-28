#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/Users/evanelko/Documents/pepseq_snakemake"
SNAKEFILE="$WORKDIR/Snakefile"
CONFIGFILE="$WORKDIR/config.yaml"

if ! command -v conda >/dev/null 2>&1; then
  echo "Error: conda not found in PATH." >&2
  exit 1
fi

if ! command -v snakemake >/dev/null 2>&1; then
  echo "Error: snakemake not found in PATH." >&2
  exit 1
fi

# Force conda environment resolution/creation to osx-64.
export CONDA_SUBDIR=osx-64

echo "Using CONDA_SUBDIR=${CONDA_SUBDIR}"
echo "Running Snakemake in: ${WORKDIR}"

cd "$WORKDIR"

CMD=(
  snakemake
  --use-conda
  -s "$SNAKEFILE"
  --configfile "$CONFIGFILE"
)

"${CMD[@]}" "$@"
