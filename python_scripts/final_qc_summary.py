#!/usr/bin/env python3
"""Build a final per-sample QC summary table across pipeline stages."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create final per-sample QC summary TSV.")
    parser.add_argument("rc_totals_tsv", type=Path, help="TSV with columns: Sample, Total Read Counts")
    parser.add_argument("rc_pass_tsv", type=Path, help="Passing sample list from RC filter")
    parser.add_argument("zero_counts_tsv", type=Path, help="TSV with columns: Sample, Zero Count")
    parser.add_argument("zero_pass_tsv", type=Path, help="Passing sample list from zero-count filter")
    parser.add_argument(
        "corr_pairs_all_tsv",
        type=Path,
        help="TSV with all pre-cutoff pair correlations (column-sum normalized stage)",
    )
    parser.add_argument("corr_pass_tsv", type=Path, help="Passing sample list from first correlation filter")
    parser.add_argument(
        "--zscore-pairs-all-tsv",
        type=Path,
        default=None,
        help="Optional TSV with all pre-cutoff pair correlations from z-score correlation QC.",
    )
    parser.add_argument(
        "--zscore-pass-tsv",
        type=Path,
        default=None,
        help="Optional passing sample list from z-score correlation QC.",
    )
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output summary TSV.")
    return parser.parse_args()


def read_pass_list(path: Path) -> set[str]:
    """Read one-sample-per-line pass-list files used by subjoin stages."""
    samples: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            sample = line.strip()
            if sample:
                samples.add(sample)
    return samples


def sample_max_correlation(path: Path) -> dict[str, float]:
    """Collapse pairwise rows into one max correlation value per sample."""
    if not path.exists() or path.stat().st_size == 0:
        return {}

    df = pd.read_csv(path, sep="\t")
    required = {"Pair1", "Pair2", "Pearson Correlation"}
    if not required.issubset(df.columns):
        return {}

    max_corr: dict[str, float] = {}
    for _, row in df.iterrows():
        corr = row["Pearson Correlation"]
        if pd.isna(corr):
            continue
        for sample_col in ("Pair1", "Pair2"):
            sample = row[sample_col]
            if sample not in max_corr or corr > max_corr[sample]:
                max_corr[sample] = float(corr)
    return max_corr


def maybe_value(condition: bool, value: float | int | None) -> float | int | str:
    """Return blank cell when a metric is not applicable downstream."""
    if not condition or value is None or pd.isna(value):
        return ""
    return value


def main() -> None:
    """Build the final wide QC table with staged blanking semantics."""
    args = parse_args()

    rc_df = pd.read_csv(args.rc_totals_tsv, sep="\t")
    if not {"Sample", "Total Read Counts"}.issubset(rc_df.columns):
        raise ValueError("RC totals file must contain columns: Sample, Total Read Counts")
    samples = rc_df["Sample"].tolist()
    rc_map = dict(zip(rc_df["Sample"], rc_df["Total Read Counts"]))

    zero_df = pd.read_csv(args.zero_counts_tsv, sep="\t")
    if not {"Sample", "Zero Count"}.issubset(zero_df.columns):
        raise ValueError("Zero-count file must contain columns: Sample, Zero Count")
    zero_map = dict(zip(zero_df["Sample"], zero_df["Zero Count"]))

    rc_pass = read_pass_list(args.rc_pass_tsv)
    zero_pass = read_pass_list(args.zero_pass_tsv)
    corr_pass = read_pass_list(args.corr_pass_tsv)

    corr_map = sample_max_correlation(args.corr_pairs_all_tsv)

    use_zscore = (
        args.zscore_pairs_all_tsv is not None
        and args.zscore_pass_tsv is not None
        and args.zscore_pairs_all_tsv.exists()
        and args.zscore_pass_tsv.exists()
    )
    zscore_pass: set[str] = set()
    zscore_map: dict[str, float] = {}
    if use_zscore:
        zscore_pass = read_pass_list(args.zscore_pass_tsv)
        zscore_map = sample_max_correlation(args.zscore_pairs_all_tsv)

    rows = []
    for sample in samples:
        rc_ok = sample in rc_pass
        zero_ok = sample in zero_pass
        corr_ok = sample in corr_pass
        zscore_ok = sample in zscore_pass if use_zscore else False

        rows.append(
            {
                "Sample": sample,
                "Read Counts": maybe_value(True, rc_map.get(sample)),
                "Zero Counts": maybe_value(rc_ok, zero_map.get(sample)),
                "Column Sum Normalized Correlation": maybe_value(zero_ok, corr_map.get(sample)),
                "Z Score Correlation": maybe_value(corr_ok and zscore_ok, zscore_map.get(sample)),
            }
        )

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"Output: {args.output}")
    print(f"Samples summarized: {len(out_df)}")


if __name__ == "__main__":
    main()
