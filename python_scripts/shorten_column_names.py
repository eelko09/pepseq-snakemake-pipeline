#!/usr/bin/env python3
"""Shorten TSV column names to the first two underscore-delimited parts.

Reads a TSV with 'Sequence name' as index column and rewrites sample column
headers by keeping only the first two underscore-separated tokens.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

SEQUENCE_COL = "Sequence name"


def parse_args() -> argparse.Namespace:
    # CLI arguments for column shortening
    parser = argparse.ArgumentParser(
        description=(
            "Shorten column names by keeping only the first two underscore-separated parts."
        )
    )
    parser.add_argument("input_tsv", type=Path, help="Input TSV path.")
    parser.add_argument("output_tsv", type=Path, help="Output TSV path.")
    parser.add_argument(
        "--map-output",
        type=Path,
        required=True,
        help="Output TSV mapping file (original name -> shortened name). Use /dev/null to skip.",
    )
    parser.add_argument(
        "--drop-prefix",
        action="append",
        default=[],
        help=(
            "Drop columns whose name starts with this prefix. "
            "May be provided multiple times."
        ),
    )
    parser.add_argument(
        "--drop-nan",
        action="store_true",
        help="Drop columns that contain any NaN values.",
    )
    return parser.parse_args()


def shorten_name(name: str) -> str:
    # Keep Sblk controls intact; otherwise use first two underscore tokens
    text = str(name)
    if text.startswith("Sblk"):
        return text
    parts = text.split("_")
    if len(parts) <= 2:
        return text
    return "_".join(parts[:2])


def main() -> None:
    # Load, shorten, and emit mapping
    args = parse_args()

    df = pd.read_csv(args.input_tsv, sep="\t", index_col=SEQUENCE_COL)
    drop_cols = []
    if args.drop_prefix:
        drop_cols.extend(
            col for col in df.columns if any(str(col).startswith(p) for p in args.drop_prefix)
        )
    if args.drop_nan:
        nan_cols = df.columns[df.isna().any()].tolist()
        drop_cols.extend(nan_cols)
    if drop_cols:
        drop_cols = sorted(set(col for col in drop_cols if col != SEQUENCE_COL))
        print("Dropped columns:")
        for col in drop_cols:
            print(col)
        if drop_cols:
            df = df.drop(columns=drop_cols)

    new_columns = []
    mapping_rows = []
    for col in df.columns:
        new_name = shorten_name(col)
        new_columns.append(new_name)
        mapping_rows.append((col, new_name))

    if len(new_columns) != len(set(new_columns)):
        duplicates = pd.Series(new_columns).duplicated(keep=False)
        dup_names = sorted(set(pd.Series(new_columns)[duplicates]))
        raise ValueError(f"Non-unique shortened names detected: {', '.join(dup_names)}")

    df.columns = new_columns
    df.to_csv(args.output_tsv, sep="\t", index=True)
    if str(args.map_output) != "/dev/null":
        mapping_df = pd.DataFrame(mapping_rows, columns=["Original Name", "Shortened Name"])
        mapping_df.to_csv(args.map_output, sep="\t", index=False)

    print(f"Input: {args.input_tsv}")
    print(f"Output: {args.output_tsv}")
    if str(args.map_output) != "/dev/null":
        print(f"Mapping: {args.map_output}")


if __name__ == "__main__":
    main()
