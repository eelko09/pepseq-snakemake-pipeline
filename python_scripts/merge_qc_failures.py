#!/usr/bin/env python3
"""Merge per-step QC failure logs into one TSV."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    # CLI arguments for merging per-step failure logs
    parser = argparse.ArgumentParser(description="Merge QC failure TSV files.")
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="Input per-step failure TSV files with columns: Sample, Failed QC Step, Failure Value",
    )
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output merged TSV")
    return parser.parse_args()


def main() -> None:
    # Read all QC failure files and combine into a single table
    args = parse_args()

    frames = []
    for path in args.inputs:
        df = pd.read_csv(path, sep="\t")
        expected = {"Sample", "Failed QC Step", "Failure Value"}
        if not expected.issubset(df.columns):
            raise ValueError(
                f"{path} missing required columns: Sample, Failed QC Step, Failure Value"
            )
        frames.append(df[["Sample", "Failed QC Step", "Failure Value"]])

    if frames:
        merged = pd.concat(frames, ignore_index=True).drop_duplicates()
    else:
        merged = pd.DataFrame(columns=["Sample", "Failed QC Step", "Failure Value"])

    merged.to_csv(args.output, sep="\t", index=False)
    print(f"Merged failures: {len(merged)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
