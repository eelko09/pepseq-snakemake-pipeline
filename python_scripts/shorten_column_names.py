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


def stream_rename_tsv(input_tsv: Path, output_tsv: Path, rename_map: dict[str, str]) -> None:
    """Rewrite a TSV in chunks while renaming columns to avoid full-table memory use."""
    first_chunk = True
    for chunk in pd.read_csv(input_tsv, sep="\t", chunksize=50000):
        chunk = chunk.rename(columns=rename_map)
        chunk.to_csv(
            output_tsv,
            sep="\t",
            index=False,
            mode="w" if first_chunk else "a",
            header=first_chunk,
        )
        first_chunk = False


def main() -> None:
    # Load, shorten, and emit mapping
    args = parse_args()

    header_df = pd.read_csv(args.input_tsv, sep="\t", nrows=0)
    if SEQUENCE_COL not in header_df.columns:
        raise ValueError(f"Input file must contain '{SEQUENCE_COL}' column.")
    sample_cols = [c for c in header_df.columns if c != SEQUENCE_COL]

    new_columns = []
    mapping_rows = []
    for col in sample_cols:
        new_name = shorten_name(col)
        new_columns.append(new_name)
        mapping_rows.append((col, new_name))

    if len(new_columns) != len(set(new_columns)):
        duplicates = pd.Series(new_columns).duplicated(keep=False)
        dup_names = sorted(set(pd.Series(new_columns)[duplicates]))
        raise ValueError(f"Non-unique shortened names detected: {', '.join(dup_names)}")

    rename_map = {orig: new for orig, new in mapping_rows}

    # Fast path used by the pipeline rule: stream rewrite with renamed headers.
    if not args.drop_prefix and not args.drop_nan:
        stream_rename_tsv(args.input_tsv, args.output_tsv, rename_map)
        if str(args.map_output) != "/dev/null":
            mapping_df = pd.DataFrame(mapping_rows, columns=["Original Name", "Shortened Name"])
            mapping_df.to_csv(args.map_output, sep="\t", index=False)
        print(f"Input: {args.input_tsv}")
        print(f"Output: {args.output_tsv}")
        if str(args.map_output) != "/dev/null":
            print(f"Mapping: {args.map_output}")
        return

    # Fallback path for optional column-dropping features.
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
