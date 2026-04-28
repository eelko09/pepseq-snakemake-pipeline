#!/usr/bin/env python3
"""Remove control sample columns from a read-count TSV.

Controls are inferred from a control matrix TSV's sample columns.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

SEQUENCE_COL = "Sequence name"


def parse_args() -> argparse.Namespace:
    # CLI arguments for control removal + name restoration
    parser = argparse.ArgumentParser(description="Remove control sample columns from raw matrix.")
    parser.add_argument("--raw", required=True, type=Path, help="Input raw count TSV")
    parser.add_argument("--controls", required=True, type=Path, help="Control matrix TSV")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV path")
    parser.add_argument("--name-map", required=True, type=Path, help="Column name map TSV")
    return parser.parse_args()


def get_sample_columns(tsv: Path) -> list[str]:
    # Extract sample columns from header, ignoring the sequence column
    header = pd.read_csv(tsv, sep="\t", nrows=0).columns.tolist()
    if not header:
        raise ValueError(f"No columns found in {tsv}")
    if SEQUENCE_COL in header:
        return [c for c in header if c != SEQUENCE_COL]
    return header[1:]


def main() -> None:
    # Load raw data, drop controls, restore full names, and write output
    args = parse_args()

    raw_df = pd.read_csv(args.raw, sep="\t")
    if SEQUENCE_COL not in raw_df.columns:
        raw_df = raw_df.rename(columns={raw_df.columns[0]: SEQUENCE_COL})

    name_map_df = pd.read_csv(args.name_map, sep="\t")
    if not {"Original Name", "Shortened Name"}.issubset(name_map_df.columns):
        raise ValueError("Name map file must contain 'Original Name' and 'Shortened Name' columns.")
    shorten_to_full = dict(zip(name_map_df["Shortened Name"], name_map_df["Original Name"]))
    raw_samples = [c for c in raw_df.columns if c != SEQUENCE_COL]
    control_samples = set()
    if str(args.controls) != "/dev/null":
        control_samples = set(get_sample_columns(args.controls))

    keep_samples = [c for c in raw_samples if c not in control_samples]
    output_df = raw_df[[SEQUENCE_COL] + keep_samples]
    # Restore full-length column names where possible.
    restored_cols = [SEQUENCE_COL]
    for col in keep_samples:
        restored_cols.append(shorten_to_full.get(col, col))
    output_df.columns = restored_cols
    output_df.to_csv(args.output, sep="\t", index=False)

    print(f"Raw samples: {len(raw_samples)}")
    print(f"Control samples to remove: {len(control_samples)}")
    print(f"Samples kept: {len(keep_samples)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
