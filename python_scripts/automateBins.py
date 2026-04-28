#!/usr/bin/env python3
"""Build control bins from superblock samples.

Pipeline steps:
1) Extract Sblk sample names from the input matrix header.
2) Subset input to Sblk-only columns.
3) Optionally filter Sblk columns by total read-count sum.
4) Normalize Sblk data with `pepsirf norm`.
5) Compute pairwise Pearson correlations and per-sample summed correlations.
6) Select top-correlated Sblks.
7) Optionally generate pairwise scatterplots for selected Sblks.
8) Subset to selected controls and build bins with `pepsirf bin`.
"""

import argparse
import logging
import subprocess
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr


LOGGER = logging.getLogger(__name__)
COMMAND_LOG_PATH = Path("automateBins.log")


class OnlyReadcountSummaryFilter(logging.Filter):
    """Allow only the read-count summary line to be emitted to the console."""

    def filter(self, record: logging.LogRecord) -> bool:
        return record.getMessage().startswith("Read count sum filter kept ")


def setup_logging() -> None:
    """Log everything to automateBins.log, but only one summary line to terminal."""
    formatter = logging.Formatter("%(levelname)s: %(message)s")
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    root_logger.handlers.clear()

    file_handler = logging.FileHandler(COMMAND_LOG_PATH, mode="w", encoding="utf-8")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    console_handler.addFilter(OnlyReadcountSummaryFilter())
    root_logger.addHandler(console_handler)


def main() -> None:
    """Run the end-to-end Sblk pipeline: extract, filter, normalize, correlate, and bin."""
    # CLI arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", help="Raw counts TSV file", required=True)
    parser.add_argument(
        "-s",
        "--sblkID",
        help="Prefix associated with all superblock samples",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Base prefix for generated output files (no extension needed).",
        required=True,
    )
    parser.add_argument(
        "--top-fraction",
        type=float,
        default=0.20,
        help="Fraction of Sblks (by summed correlation) to include for bins.",
    )
    parser.add_argument(
        "--min-readcount-sum",
        type=int,
        default=0,
        help="Minimum allowed sum of counts for a sample column to be kept.",
    )
    parser.add_argument(
        "--no-scatterplots",
        action="store_true",
        help="Disable generation of pairwise Sblk scatterplots.",
    )
    parser.add_argument(
        "--max-pairs",
        type=int,
        default=None,
        help="Maximum number of Sblk pairs to plot. Default plots all pairs.",
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
    args = parser.parse_args()

    # Logging and basic validation
    setup_logging()

    input_path = Path(args.input)
    if not input_path.exists():
        raise SystemExit(f"Input file not found: {input_path}")

    if not 0 < args.top_fraction <= 1:
        raise SystemExit("--top-fraction must be in the range (0, 1].")
    if args.min_readcount_sum < 0:
        raise SystemExit("--min-readcount-sum must be >= 0.")

    output_prefix = args.outfile
    bins_dir = Path("bins")
    scatterplot_dir = bins_dir / "Sblk_scatterplots"
    make_scatterplots = not args.no_scatterplots

    bins_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Directory '%s' is ready.", bins_dir)

    # 1) Identify superblock (Sblk) sample columns from the original input table.
    sblk_names_path = Path("sblk_names.txt")
    extract_sblk_sample_names(input_path, args.sblkID, sblk_names_path)

    # 2) Extract only Sblk samples from raw counts.
    sblk_raw_counts_file = f"{output_prefix}_sblkOnly.tsv"
    sblk_raw_counts_path = bins_dir / sblk_raw_counts_file
    run_linux_command(
        "pepsirf",
        "subjoin",
        "-i",
        f"{input_path},{sblk_names_path}",
        "-o",
        str(sblk_raw_counts_path),
    )

    # 3) Filter extracted Sblk samples by total read-count sum threshold.
    sblk_raw_counts_filtered_path = bins_dir / f"{output_prefix}_sblkOnly_filtered.tsv"
    filter_samples_by_total_counts(
        sblk_raw_counts_path,
        sblk_raw_counts_filtered_path,
        min_readcount_sum=args.min_readcount_sum,
    )

    # 4) Column-sum normalize filtered Sblk counts.
    sblk_csnorm_file = f"{output_prefix}_sblkOnly_CS.tsv"
    run_linux_command(
        "pepsirf",
        "norm",
        "-p",
        str(sblk_raw_counts_filtered_path),
        "-o",
        str(bins_dir / sblk_csnorm_file),
    )

    # 5) Compute all pairwise Sblk Pearson correlations for ranking.
    sblk_cs_df = pd.read_csv(
        bins_dir / sblk_csnorm_file, sep="\t", index_col="Sequence name"
    )
    pearson_correlations_with_sums(
        sblk_cs_df,
        bins_dir / f"{output_prefix}_Sblk_correlations.tsv",
        bins_dir / f"{output_prefix}_Sblk_correlations_summed.tsv",
    )

    # 6) Select top correlated Sblks to use as controls for binning.
    sblk_corrs_path = bins_dir / f"{output_prefix}_Sblk_correlations_summed.tsv"
    sblk_corrs_df = pd.read_csv(sblk_corrs_path, sep="\t")
    if sblk_corrs_df.empty:
        raise SystemExit(
            "No valid Sblk correlations were produced; cannot continue to binning."
        )

    rows_to_use = max(1, int(len(sblk_corrs_df) * args.top_fraction))
    sblk_corrs_df_sorted = sblk_corrs_df.sort_values(
        by="SumOfCorrelations", ascending=False
    )
    cutoff_value = sblk_corrs_df_sorted.iloc[rows_to_use - 1]["SumOfCorrelations"]
    if cutoff_value > (0.1 * len(sblk_corrs_df_sorted)):
        sblk_corrs_df_to_use = sblk_corrs_df_sorted.head(rows_to_use)
        sblk_names_for_bins_path = Path("Sblk_names_forBins.txt")
        sblk_corrs_df_to_use["Sample"].to_csv(
            sblk_names_for_bins_path, sep="\t", index=False
        )
        selected_samples = sblk_corrs_df_to_use["Sample"].tolist()
    else:
        LOGGER.warning(
            "Low correlation of Sblks; continuing with top %.0f%% anyway.",
            args.top_fraction * 100,
        )
        sblk_names_for_bins_path = Path("Sblk_names_forBins.txt")
        sblk_corrs_df_sorted.head(rows_to_use)["Sample"].to_csv(
            sblk_names_for_bins_path, sep="\t", index=False
        )
        selected_samples = sblk_corrs_df_sorted.head(rows_to_use)["Sample"].tolist()

    # 7) Create scatterplots only for the final selected Sblk samples.
    if make_scatterplots:
        scatterplot_dir.mkdir(parents=True, exist_ok=True)
        LOGGER.info("Directory '%s' is ready.", scatterplot_dir)
        plot_pairwise_scatterplots(
            sblk_cs_df[selected_samples],
            scatterplot_dir=scatterplot_dir,
            plot_format=args.plot_format,
            plot_dpi=args.plot_dpi,
            max_pairs=args.max_pairs,
        )

    # 8) Subset normalized Sblks to selected controls, then build bins.
    sblk_cs_bins_file = f"{output_prefix}_sblkOnly_CS_forBins.tsv"
    run_linux_command(
        "pepsirf",
        "subjoin",
        "-i",
        f"{bins_dir / sblk_csnorm_file},{sblk_names_for_bins_path}",
        "-o",
        str(bins_dir / sblk_cs_bins_file),
    )

    bins_file = f"{output_prefix}_b300r1_bins.tsv"
    run_linux_command(
        "pepsirf",
        "bin",
        "-s",
        str(bins_dir / sblk_cs_bins_file),
        "-b",
        "300",
        "-r",
        "1",
        "-o",
        str(bins_dir / bins_file),
    )

    LOGGER.info("Binning complete.")


def extract_sblk_sample_names(tsv_path: Path, prefix: str, output_path: Path) -> None:
    """Write all column names starting with prefix (e.g., 'Sblk') to a one-per-line text file."""
    with open(tsv_path, "r", encoding="utf-8") as infile:
        header_line = infile.readline().strip()
        columns = header_line.split("\t")
        matching_columns = [col for col in columns if col.startswith(prefix)]

    if not matching_columns:
        raise SystemExit(
            f"No sample columns found with prefix '{prefix}' in file: {tsv_path}"
        )

    with open(output_path, "w", encoding="utf-8") as outfile:
        for col in matching_columns:
            outfile.write(col + "\n")


def filter_samples_by_total_counts(
    input_tsv_path: Path, output_tsv_path: Path, min_readcount_sum: int
) -> None:
    """Keep only sample columns whose column-sum is >= min_readcount_sum; preserve first ID column."""
    df = pd.read_csv(input_tsv_path, sep="\t")
    if df.shape[1] < 2:
        raise SystemExit(
            f"Input file must contain an identifier column and at least one sample column: {input_tsv_path}"
        )

    id_column = df.columns[0]
    sample_columns = list(df.columns[1:])
    numeric_samples = df[sample_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    sample_sums = numeric_samples.sum(axis=0)
    kept_columns = [
        col for col in sample_columns if sample_sums[col] >= min_readcount_sum
    ]

    if not kept_columns:
        raise SystemExit(
            f"No sample columns meet --min-readcount-sum {min_readcount_sum} in file: {input_tsv_path}"
        )

    filtered_df = pd.concat([df[[id_column]], df[kept_columns]], axis=1)
    filtered_df.to_csv(output_tsv_path, sep="\t", index=False)
    LOGGER.info(
        "Read count sum filter kept %d/%d sample columns (threshold=%s). Refined file: %s",
        len(kept_columns),
        len(sample_columns),
        min_readcount_sum,
        output_tsv_path,
    )


def pearson_correlations_with_sums(
    df: pd.DataFrame,
    pairwise_output_path: Path,
    sum_output_path: Path,
) -> None:
    """
    Compute all pairwise Pearson correlations between sample columns.
    Save pairwise and summed-correlation TSV outputs for sample ranking.
    """
    results = []
    correlation_sums = defaultdict(float)

    numeric_df = df.apply(pd.to_numeric, errors="coerce")

    for col1, col2 in combinations(numeric_df.columns, 2):
        aligned = numeric_df[[col1, col2]].dropna()
        if len(aligned) < 2:
            continue

        # Pearson is undefined for constant vectors.
        if aligned[col1].nunique() < 2 or aligned[col2].nunique() < 2:
            continue

        r, _ = pearsonr(aligned[col1], aligned[col2])
        if pd.isna(r):
            continue

        results.append((col1, col2, r))
        correlation_sums[col1] += r
        correlation_sums[col2] += r

    results_sorted = sorted(results, key=lambda row: row[2], reverse=True)
    correlation_sums_sorted = sorted(
        correlation_sums.items(), key=lambda item: item[1], reverse=True
    )

    with open(pairwise_output_path, "w", encoding="utf-8") as f:
        f.write("Sample1\tSample2\tPearsonCorrelation\n")
        for col1, col2, r in results_sorted:
            f.write(f"{col1}\t{col2}\t{r:.6f}\n")

    with open(sum_output_path, "w", encoding="utf-8") as f:
        f.write("Sample\tSumOfCorrelations\n")
        for col, total_r in correlation_sums_sorted:
            f.write(f"{col}\t{total_r:.6f}\n")

    LOGGER.info("Computed %d pairwise correlations.", len(results))


def plot_pairwise_scatterplots(
    df: pd.DataFrame,
    scatterplot_dir: Path = Path("bins/Sblk_scatterplots"),
    plot_format: str = "png",
    plot_dpi: int = 300,
    max_pairs: int | None = None,
) -> None:
    """Create x=y scatterplots for all pairwise sample comparisons in df."""
    plotted_pairs = 0
    numeric_df = df.apply(pd.to_numeric, errors="coerce")

    for col1, col2 in combinations(numeric_df.columns, 2):
        if max_pairs is not None and plotted_pairs >= max_pairs:
            break

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
        plotted_pairs += 1

    LOGGER.info("Generated %d scatterplots for selected Sblk samples.", plotted_pairs)


def sanitize_filename(name: str) -> str:
    """Convert arbitrary sample names into filesystem-safe plot filename fragments."""
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.")
    return "".join(char if char in allowed else "_" for char in str(name))


def run_linux_command(*inputs: object) -> None:
    """Run an external command and stop execution with a clear error message on failure."""
    if not inputs:
        raise ValueError("At least one command input must be provided.")

    cmd = [str(part) for part in inputs]
    try:
        with open(COMMAND_LOG_PATH, "a", encoding="utf-8") as command_log:
            command_log.write(f"Running: {' '.join(cmd)}\n")
            subprocess.run(cmd, check=True, stdout=command_log, stderr=command_log)
    except FileNotFoundError as exc:
        raise SystemExit(f"Required executable not found: {cmd[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise SystemExit(
            f"Command failed with exit code {exc.returncode}: {' '.join(cmd)}"
        ) from exc


if __name__ == "__main__":
    main()
