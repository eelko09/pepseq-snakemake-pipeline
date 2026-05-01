#!/usr/bin/env python3
"""Compute pairwise Pearson correlations from a pairs file.

Inputs:
1) Pairs TSV with columns: Pair1, Pair2, Comparison
2) Data TSV with index column: Sequence name
   and sample columns matching Pair1/Pair2 names.

Outputs:
1) Filtered pairwise correlation TSV with columns:
   Pair1, Pair2, Pearson Correlation, Comparison
2) Passing sample-name TSV (one column, no header) for downstream subjoin.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


REQUIRED_PAIR_COLS = {"Pair1", "Pair2", "Comparison"}
OUTPUT_COLS = [
    "Pair1",
    "Pair2",
    "Pearson Correlation",
    "Pearson Correlation (log2(Z+8)-3)",
    "Comparison",
]
SEQUENCE_COL = "Sequence name"
STEP_NAME = "CORRELATION_FILTER"
DEFAULT_ALWAYS_KEEP_PREFIX = "Sblk"


def parse_args() -> argparse.Namespace:
    # CLI arguments for pairwise correlation filtering
    parser = argparse.ArgumentParser(
        description="Compute Pearson correlation for each Pair1/Pair2 row."
    )
    parser.add_argument("pairs_tsv", type=Path, help="Path to pairs TSV file.")
    parser.add_argument(
        "data_tsv",
        type=Path,
        help="Path to data TSV file with 'Sequence name' as index.",
    )
    parser.add_argument(
        "cutoff",
        type=float,
        help="Pearson correlation cutoff. Keeps rows where correlation >= cutoff.",
    )
    parser.add_argument(
        "--log-pearson-cutoff",
        type=float,
        default=None,
        help=(
            "Optional cutoff for Pearson correlation computed from log-transformed "
            "values: log2(value + 8) - 3. When set, samples must also pass this "
            "cutoff to be retained."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help=(
            "Output TSV path. Defaults to <pairs_file_stem>_pearson.tsv in the "
            "pairs file directory."
        ),
    )
    parser.add_argument(
        "--passing-samples-output",
        type=Path,
        default=None,
        help=(
            "Output path for passing sample names (single column, no header). "
            "Defaults to <data_file_stem>_correlation_pass_samples.tsv."
        ),
    )
    parser.add_argument(
        "--failed-samples-output",
        type=Path,
        default=None,
        help="Optional TSV output path for samples removed by this QC step.",
    )
    parser.add_argument(
        "--all-correlations-output",
        type=Path,
        default=None,
        help=(
            "Optional TSV output path for all computed pair correlations "
            "(before cutoff filtering)."
        ),
    )
    parser.add_argument(
        "--always-keep-prefix",
        type=str,
        default=DEFAULT_ALWAYS_KEEP_PREFIX,
        help=(
            "Sample-name prefix to always keep in filtered data output, even if not "
            "retained by pairwise correlation filtering."
        ),
    )
    parser.add_argument(
        "--keep-samples-not-in-pairs",
        action="store_true",
        help=(
            "Also retain input samples that do not appear in Pair1/Pair2 anywhere in "
            "the pairs file. Samples that do appear in pairs must still pass correlation QC."
        ),
    )
    parser.add_argument(
        "--scatterplots",
        action="store_true",
        help="Generate scatterplots for kept pairs.",
    )
    parser.add_argument(
        "--scatterplot-dir",
        type=Path,
        default=Path("pairwise_scatterplots"),
        help="Directory for scatterplot outputs.",
    )
    parser.add_argument(
        "--plot-format",
        choices=("png", "pdf", "svg"),
        default="png",
        help="Output format for scatterplots.",
    )
    parser.add_argument(
        "--plot-dpi",
        type=int,
        default=300,
        help="DPI for raster scatterplots.",
    )
    parser.add_argument(
        "--max-pairs",
        type=int,
        default=None,
        help="Maximum number of pairs to plot. Default plots all kept pairs.",
    )
    return parser.parse_args()


def default_output_path(pairs_path: Path) -> Path:
    # Place results alongside the pairs file by default
    suffix = pairs_path.suffix if pairs_path.suffix else ".tsv"
    return pairs_path.with_name(f"{pairs_path.stem}_pearson{suffix}")


def default_all_correlations_output_path(pairs_path: Path) -> Path:
    # Place pre-cutoff correlations alongside the pairs file by default
    suffix = pairs_path.suffix if pairs_path.suffix else ".tsv"
    return pairs_path.with_name(f"{pairs_path.stem}_pearson_all{suffix}")


def default_passing_samples_output_path(data_path: Path) -> Path:
    # Place pass-list output alongside the input data by default
    suffix = data_path.suffix if data_path.suffix else ".tsv"
    return data_path.with_name(f"{data_path.stem}_correlation_pass_samples{suffix}")


def validate_pairs_columns(pairs_df: pd.DataFrame) -> None:
    # Fail fast if the pairs file is not in the expected format
    missing = REQUIRED_PAIR_COLS.difference(pairs_df.columns)
    if missing:
        missing_list = ", ".join(sorted(missing))
        raise ValueError(f"Pairs file is missing required columns: {missing_list}")


def transform_for_log_correlation(data_df: pd.DataFrame) -> pd.DataFrame:
    # Apply log transform used for secondary correlation metric.
    shifted = data_df + 8.0
    shifted = shifted.where(shifted > 0.0)
    return np.log2(shifted) - 3.0


def compute_correlations(
    pairs_df: pd.DataFrame, data_df: pd.DataFrame, log_data_df: pd.DataFrame
) -> pd.DataFrame:
    # Compute Pearson r for each pair, skipping missing columns
    results = []
    skipped_missing = 0
    corr_cache = {}
    log_corr_cache = {}

    for i, row in pairs_df.iterrows():
        pair1 = row["Pair1"]
        pair2 = row["Pair2"]
        comparison = row["Comparison"]

        if pair1 not in data_df.columns:
            skipped_missing += 1
            print(
                f"Warning: row {i} skipped; Pair1 column '{pair1}' not found in data TSV",
                file=sys.stderr,
            )
            continue
        if pair2 not in data_df.columns:
            skipped_missing += 1
            print(
                f"Warning: row {i} skipped; Pair2 column '{pair2}' not found in data TSV",
                file=sys.stderr,
            )
            continue

        # Cache symmetric correlations to avoid recomputing duplicate pairs
        cache_key = (pair1, pair2) if pair1 <= pair2 else (pair2, pair1)
        pearson_r = corr_cache.get(cache_key)
        if pearson_r is None:
            # Pearson correlation uses pairwise complete observations.
            pearson_r = data_df[pair1].corr(data_df[pair2], method="pearson")
            corr_cache[cache_key] = pearson_r
        log_pearson_r = log_corr_cache.get(cache_key)
        if log_pearson_r is None:
            log_pearson_r = log_data_df[pair1].corr(log_data_df[pair2], method="pearson")
            log_corr_cache[cache_key] = log_pearson_r

        results.append(
            {
                "Pair1": pair1,
                "Pair2": pair2,
                "Pearson Correlation": pearson_r,
                "Pearson Correlation (log2(Z+8)-3)": log_pearson_r,
                "Comparison": comparison,
            }
        )

    out_df = pd.DataFrame(results, columns=OUTPUT_COLS)
    if skipped_missing:
        print(f"Skipped pairs with missing columns: {skipped_missing}", file=sys.stderr)
    return out_df


def build_keep_columns(
    filtered_pairs_df: pd.DataFrame,
    always_keep_samples: list[str],
    input_sample_order: list[str],
    extra_keep_samples: list[str] | None = None,
) -> list[str]:
    # Preserve original column order while keeping selected samples
    keep_set = set(always_keep_samples)
    if extra_keep_samples is not None:
        keep_set.update(extra_keep_samples)

    for _, row in filtered_pairs_df.iterrows():
        for col in (row["Pair1"], row["Pair2"]):
            keep_set.add(col)

    return [col for col in input_sample_order if col in keep_set]


def write_failed_samples_tsv(
    path: Path, failed_samples: list[tuple[str, object]], step_name: str
) -> None:
    # Record samples removed by this step with their failure values
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "Failed QC Step", "Failure Value"])
        for sample in failed_samples:
            writer.writerow([sample[0], step_name, sample[1]])


def write_passing_samples_tsv(path: Path, passing_samples: list[str]) -> None:
    # Write one passing sample per line for downstream pepsirf subjoin
    with path.open("w", encoding="utf-8") as handle:
        for sample in passing_samples:
            handle.write(f"{sample}\n")


def build_sample_max_correlations(
    corr_df: pd.DataFrame, corr_column: str = "Pearson Correlation"
) -> dict[str, float]:
    # Reduce multiple pair correlations to a single max per sample
    sample_max_corr: dict[str, float] = {}
    for _, row in corr_df.iterrows():
        r = row[corr_column]
        if pd.isna(r):
            continue
        p1 = row["Pair1"]
        p2 = row["Pair2"]
        if p1 not in sample_max_corr or r > sample_max_corr[p1]:
            sample_max_corr[p1] = r
        if p2 not in sample_max_corr or r > sample_max_corr[p2]:
            sample_max_corr[p2] = r
    return sample_max_corr


def sanitize_filename(name: str) -> str:
    # Make sample names safe for filenames
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.")
    return "".join(char if char in allowed else "_" for char in str(name))


def plot_pairwise_scatterplots(
    data_df: pd.DataFrame,
    pairs_df: pd.DataFrame,
    scatterplot_dir: Path,
    plot_format: str,
    plot_dpi: int,
    max_pairs: int | None,
) -> int:
    # Optional scatterplot generation for kept pairs
    scatterplot_dir.mkdir(parents=True, exist_ok=True)
    plotted = 0

    numeric_df = data_df.apply(pd.to_numeric, errors="coerce")
    for _, row in pairs_df.iterrows():
        if max_pairs is not None and plotted >= max_pairs:
            break
        col1 = row["Pair1"]
        col2 = row["Pair2"]
        if col1 not in numeric_df.columns or col2 not in numeric_df.columns:
            continue
        aligned = numeric_df[[col1, col2]].dropna()
        if len(aligned) < 2:
            continue

        plt.figure(figsize=(6, 6))
        sns.scatterplot(data=aligned, x=col1, y=col2, s=14, alpha=0.7)
        min_v = min(aligned[col1].min(), aligned[col2].min())
        max_v = max(aligned[col1].max(), aligned[col2].max())
        plt.plot(
            [min_v, max_v],
            [min_v, max_v],
            linestyle="--",
            linewidth=1,
            color="black",
            label="x=y",
        )
        plt.xlabel(col1)
        plt.ylabel(col2)
        plt.title(f"{col1} vs {col2}")
        plt.legend(loc="best")
        plt.tight_layout()
        output_plot = (
            scatterplot_dir
            / f"{sanitize_filename(col1)}_vs_{sanitize_filename(col2)}.{plot_format}"
        )
        plt.savefig(output_plot, dpi=plot_dpi)
        plt.close()
        plotted += 1

    return plotted


def main() -> None:
    # End-to-end correlation QC driver
    args = parse_args()

    # 1) Load and validate the pairs file
    pairs_df = pd.read_csv(args.pairs_tsv, sep="\t")
    validate_pairs_columns(pairs_df)

    # 2) Read sample columns and decide which columns are needed
    all_columns_df = pd.read_csv(args.data_tsv, sep="\t", nrows=0)
    if SEQUENCE_COL not in all_columns_df.columns:
        raise ValueError(f"Data file must contain '{SEQUENCE_COL}' column.")
    input_sample_columns = [c for c in all_columns_df.columns if c != SEQUENCE_COL]
    always_keep_samples = [
        c for c in input_sample_columns if c.startswith(args.always_keep_prefix)
    ]
    samples_in_pairs = set(pairs_df["Pair1"]).union(set(pairs_df["Pair2"]))
    unpaired_samples = [c for c in input_sample_columns if c not in samples_in_pairs]
    extra_keep_samples = unpaired_samples if args.keep_samples_not_in_pairs else []

    needed_columns = (
        set(pairs_df["Pair1"]).union(set(pairs_df["Pair2"])).union(set(always_keep_samples))
    )
    usecols = lambda c: c == SEQUENCE_COL or c in needed_columns
    # 3) Load only required columns for memory efficiency
    data_df = pd.read_csv(args.data_tsv, sep="\t", usecols=usecols, index_col="Sequence name")
    numeric_data_df = data_df.apply(pd.to_numeric, errors="coerce")
    log_data_df = transform_for_log_correlation(numeric_data_df)

    # 4) Compute correlations for requested pairs
    output_df = compute_correlations(pairs_df, numeric_data_df, log_data_df)
    all_correlations_output_path = (
        args.all_correlations_output
        if args.all_correlations_output is not None
        else default_all_correlations_output_path(args.pairs_tsv)
    )
    output_df.to_csv(all_correlations_output_path, sep="\t", index=False)
    no_pairs_reason = (
        "No pairs found in input data matrix (no Pair1/Pair2 rows had both samples present)."
    )
    if output_df.empty:
        # No valid pair columns were present; keep only always-keep samples.
        print(f"ERROR: {no_pairs_reason}", file=sys.stderr)
        output_path = (
            args.output if args.output is not None else default_output_path(args.pairs_tsv)
        )
        passing_samples_output_path = (
            args.passing_samples_output
            if args.passing_samples_output is not None
            else default_passing_samples_output_path(args.data_tsv)
        )
        output_df.to_csv(output_path, sep="\t", index=False)
        passing_samples = [
            c
            for c in input_sample_columns
            if c in set(always_keep_samples).union(set(extra_keep_samples))
        ]
        write_passing_samples_tsv(passing_samples_output_path, passing_samples)
        if args.failed_samples_output is not None:
            failed_samples_with_values = [
                (sample, no_pairs_reason)
                for sample in input_sample_columns
                if sample not in set(passing_samples)
            ]
            write_failed_samples_tsv(
                args.failed_samples_output, failed_samples_with_values, STEP_NAME
            )

        print(f"Pairs file: {args.pairs_tsv}")
        print(f"Data file: {args.data_tsv}")
        print(f"Cutoff: {args.cutoff}")
        print(f"Pairs processed: {len(output_df)}")
        print(f"Pairs kept: 0")
        print(f"All correlations output: {all_correlations_output_path}")
        print(f"Correlation output: {output_path}")
        print(f"Passing samples output: {passing_samples_output_path}")
        print(f"Always-kept '{args.always_keep_prefix}' samples: {len(always_keep_samples)}")
        if args.keep_samples_not_in_pairs:
            print(f"Unpaired samples auto-kept: {len(extra_keep_samples)}")
        if args.failed_samples_output is not None:
            print(f"Samples failed {STEP_NAME}: {len(failed_samples_with_values)}")
        return

    # 5) Use per-sample max correlation for cutoff when samples appear in multiple pairs.
    sample_max_corr = build_sample_max_correlations(output_df)
    sample_max_log_corr = build_sample_max_correlations(
        output_df, corr_column="Pearson Correlation (log2(Z+8)-3)"
    )

    def pair_passes(row):
        r1 = sample_max_corr.get(row["Pair1"], float("-inf"))
        r2 = sample_max_corr.get(row["Pair2"], float("-inf"))
        raw_pass = max(r1, r2) >= args.cutoff
        if args.log_pearson_cutoff is None:
            return raw_pass
        lr1 = sample_max_log_corr.get(row["Pair1"], float("-inf"))
        lr2 = sample_max_log_corr.get(row["Pair2"], float("-inf"))
        return raw_pass and (max(lr1, lr2) >= args.log_pearson_cutoff)

    # 6) Filter pairs and derive passing sample names
    pass_mask = output_df.apply(pair_passes, axis=1)
    filtered_pairs_df = output_df[pass_mask].copy()
    # Keep the union of samples present in passing pairs + always-keep samples.
    passing_samples = build_keep_columns(
        filtered_pairs_df,
        always_keep_samples,
        input_sample_columns,
        extra_keep_samples=extra_keep_samples,
    )

    output_path = args.output if args.output is not None else default_output_path(args.pairs_tsv)
    passing_samples_output_path = (
        args.passing_samples_output
        if args.passing_samples_output is not None
        else default_passing_samples_output_path(args.data_tsv)
    )
    filtered_pairs_df.to_csv(output_path, sep="\t", index=False)
    write_passing_samples_tsv(passing_samples_output_path, passing_samples)
    if args.failed_samples_output is not None:
        kept = set(passing_samples)
        failed_samples = [c for c in input_sample_columns if c not in kept]
        if args.log_pearson_cutoff is None:
            failed_value_map = build_sample_max_correlations(output_df)
        else:
            failed_value_map = build_sample_max_correlations(
                output_df, corr_column="Pearson Correlation (log2(Z+8)-3)"
            )
        failed_samples_with_values = [(s, failed_value_map.get(s, float("nan"))) for s in failed_samples]
        write_failed_samples_tsv(
            args.failed_samples_output, failed_samples_with_values, STEP_NAME
        )

    print(f"Pairs file: {args.pairs_tsv}")
    print(f"Data file: {args.data_tsv}")
    print(f"Cutoff: {args.cutoff}")
    if args.log_pearson_cutoff is not None:
        print(f"Log-transform cutoff: {args.log_pearson_cutoff}")
    print(f"Pairs processed: {len(output_df)}")
    print(f"Pairs kept: {len(filtered_pairs_df)}")
    print(f"All correlations output: {all_correlations_output_path}")
    print(f"Correlation output: {output_path}")
    print(f"Passing samples output: {passing_samples_output_path}")
    print(f"Always-kept '{args.always_keep_prefix}' samples: {len(always_keep_samples)}")
    if args.keep_samples_not_in_pairs:
        print(f"Unpaired samples auto-kept: {len(extra_keep_samples)}")
    if args.failed_samples_output is not None:
        print(f"Samples failed {STEP_NAME}: {len(failed_samples)}")

    # 7) Optional scatterplots for retained pairs
    if args.scatterplots:
        scatter_df = data_df.loc[:, passing_samples]
        plotted = plot_pairwise_scatterplots(
            scatter_df,
            filtered_pairs_df,
            scatterplot_dir=args.scatterplot_dir,
            plot_format=args.plot_format,
            plot_dpi=args.plot_dpi,
            max_pairs=args.max_pairs,
        )
        print(f"Scatterplots written: {plotted}")


if __name__ == "__main__":
    main()
