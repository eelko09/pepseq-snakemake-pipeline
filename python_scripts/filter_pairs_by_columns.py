#!/usr/bin/env python3
"""Filter a pairs TSV to only include pairs present in a data matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED_COLS = {"Pair1", "Pair2"}
SEQUENCE_COL = "Sequence name"


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for matrix-aware pair filtering."""
    parser = argparse.ArgumentParser(
        description=(
            "Subset pairs TSV to rows where Pair1 and Pair2 are both columns in the data TSV."
        )
    )
    parser.add_argument("pairs_tsv", type=Path, help="Input pairs TSV.")
    parser.add_argument("data_tsv", type=Path, help="Data TSV with sample columns.")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output filtered pairs TSV.")
    return parser.parse_args()


def main() -> None:
    """Keep only pair rows whose samples both exist in the matrix header."""
    args = parse_args()

    pairs_df = pd.read_csv(args.pairs_tsv, sep="\t")
    missing = REQUIRED_COLS.difference(pairs_df.columns)
    if missing:
        raise ValueError(f"Pairs TSV missing required columns: {', '.join(sorted(missing))}")

    header_df = pd.read_csv(args.data_tsv, sep="\t", nrows=0)
    if SEQUENCE_COL in header_df.columns:
        sample_cols = [c for c in header_df.columns if c != SEQUENCE_COL]
    else:
        sample_cols = header_df.columns.tolist()[1:]
    sample_set = set(sample_cols)

    filtered = pairs_df[
        pairs_df["Pair1"].isin(sample_set) & pairs_df["Pair2"].isin(sample_set)
    ].copy()

    filtered.to_csv(args.output, sep="\t", index=False)
    print(f"Input pairs: {len(pairs_df)}")
    print(f"Filtered pairs: {len(filtered)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
