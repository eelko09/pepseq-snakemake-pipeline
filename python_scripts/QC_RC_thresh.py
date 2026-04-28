#!/usr/bin/env python3
"""Select sample names by column-sum cutoff.

Reads a TSV file using "Sequence name" as index, sums each sample column,
and writes a no-header single-column TSV of sample names whose sums are
strictly greater than the cutoff.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd

SEQUENCE_COL = "Sequence name"
STEP_NAME = "RC_FILTER"


def parse_args() -> argparse.Namespace:
    # CLI arguments for large TSV processing
    parser = argparse.ArgumentParser(
        description=(
            "Write sample names whose column sums are strictly greater than "
            "the provided cutoff."
        )
    )
    parser.add_argument(
        "input_tsv",
        type=Path,
        help="Path to input TSV file with 'Sequence name' as index column.",
    )
    parser.add_argument(
        "cutoff",
        type=float,
        help="Column-sum cutoff. Only columns with sum > cutoff are kept.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help=(
            "Output TSV path for passing sample names (single column, no header). "
            "Defaults to <input_stem>_<cutoff>_RC_pass_samples.tsv."
        ),
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=50000,
        help="Rows per chunk for large-file processing. [50000]",
    )
    parser.add_argument(
        "--failed-samples-output",
        type=Path,
        default=None,
        help="Optional TSV output path for samples removed by this QC step.",
    )
    parser.add_argument(
        "--all-samples-output",
        type=Path,
        default=None,
        help=(
            "Optional TSV output path for total read counts of all samples. "
            "Columns: Sample, Total Read Counts."
        ),
    )
    return parser.parse_args()


def make_output_path(input_path: Path, cutoff: float) -> Path:
    # Use a stable cutoff label in the filename
    if cutoff.is_integer():
        cutoff_label = str(int(cutoff))
    else:
        cutoff_label = str(cutoff).replace(".", "p")

    return input_path.with_name(f"{input_path.stem}_{cutoff_label}_RC_pass_samples.tsv")


def compute_column_sums(input_tsv: Path, chunksize: int) -> pd.Series:
    # Stream through file in chunks to limit memory use
    sums = None
    for chunk in pd.read_csv(input_tsv, sep="\t", chunksize=chunksize):
        if SEQUENCE_COL not in chunk.columns:
            raise ValueError(f"Input file must contain '{SEQUENCE_COL}' column.")
        sample_chunk = chunk.drop(columns=[SEQUENCE_COL], errors="ignore")
        numeric_chunk = sample_chunk.apply(pd.to_numeric, errors="coerce")
        chunk_sums = numeric_chunk.sum(axis=0, skipna=True)
        sums = chunk_sums if sums is None else sums.add(chunk_sums, fill_value=0)
    if sums is None:
        raise ValueError("Input TSV is empty.")
    return sums


def write_pass_samples(path: Path, samples: list[str]) -> None:
    # One sample name per line, no header
    with path.open("w", encoding="utf-8") as handle:
        for sample in samples:
            handle.write(f"{sample}\n")


def write_failed_samples_tsv(
    path: Path, failed_samples: list[str], step_name: str, sample_values: dict[str, float]
) -> None:
    # Record samples removed by this step with their column sums
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "Failed QC Step", "Failure Value"])
        for sample in failed_samples:
            writer.writerow([sample, step_name, sample_values[sample]])


def write_all_samples_tsv(path: Path, column_sums: pd.Series) -> None:
    # Record total read counts for all samples
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "Total Read Counts"])
        for sample in column_sums.index:
            writer.writerow([sample, float(column_sums[sample])])


def main() -> None:
    # Orchestrate RC filtering
    args = parse_args()

    if args.chunksize <= 0:
        raise ValueError("--chunksize must be > 0")

    column_sums = compute_column_sums(args.input_tsv, args.chunksize)
    keep_columns = list(column_sums[column_sums > args.cutoff].index)
    keep_set = set(keep_columns)
    failed_samples = [col for col in column_sums.index if col not in keep_set]

    output_path = args.output if args.output is not None else make_output_path(args.input_tsv, args.cutoff)
    write_pass_samples(output_path, keep_columns)
    if args.failed_samples_output is not None:
        sample_values = {sample: float(column_sums[sample]) for sample in failed_samples}
        write_failed_samples_tsv(
            args.failed_samples_output, failed_samples, STEP_NAME, sample_values
        )
    if args.all_samples_output is not None:
        write_all_samples_tsv(args.all_samples_output, column_sums)

    print(f"Input: {args.input_tsv}")
    print(f"Cutoff: {args.cutoff}")
    print(f"Samples passing threshold: {len(keep_columns)} / {len(column_sums)}")
    print(f"Samples failed {STEP_NAME}: {len(failed_samples)}")
    print(f"Output: {output_path}")
    if args.all_samples_output is not None:
        print(f"All-sample read counts output: {args.all_samples_output}")


if __name__ == "__main__":
    main()
