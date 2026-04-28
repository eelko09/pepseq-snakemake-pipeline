#!/usr/bin/env python3
"""Downsample a count matrix and filter samples by zero-count distribution.

Behavior:
1) For each sample column, downsample to `--max` reads when needed.
2) Compute zero-counts per sample on the downsampled matrix.
3) Keep samples with zero-count <= mean + (`--sd-mult` * population std-dev).
4) Emit pass-list for downstream filtering plus optional QC artifacts.
"""

import csv
import optparse
import random
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd

try:
    import matrixtools as mt  # Available at https://github.com/jtladner/Modules
except ModuleNotFoundError:
    mt = None

# This script downsamples a full raw count PepSeq matrix and then:
# 1) Counts zeros in each sample column of the downsampled matrix
# 2) Calculates the zero-count distribution (mean and standard deviation)
# 3) Outputs a no-header single-column file of sample names that pass:
#    mean +/- (sd_mult * standard deviation)
STEP_NAME = "ZERO_COUNT_FILTER"


def parse_counts_fallback(path: str, delim: str = "\t") -> dict[str, dict[str, int]]:
    """Parse count matrix using pandas when matrixtools is unavailable.

    Returns a mapping of sample -> {peptide -> count} to match matrixtools.parseCounts.
    """
    df = pd.read_csv(path, sep=delim)
    if df.empty:
        return {}
    sequence_col = "Sequence name"
    if sequence_col not in df.columns:
        sequence_col = df.columns[0]
    sample_cols = [c for c in df.columns if c != sequence_col]
    counts_df = df.set_index(sequence_col)[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    out: dict[str, dict[str, int]] = {}
    for sample in sample_cols:
        out[sample] = {
            str(peptide): int(value) for peptide, value in counts_df[sample].items()
        }
    return out

def main() -> None:
    # CLI options
    usage = '%prog [options]'
    p = optparse.OptionParser(usage=usage)

    p.add_option('-d', '--data',  help='Matrix containing reads counts to downsample. [None, REQ]')
    p.add_option('-o', '--out',  help='Output file containing passing sample names, one per line and no header. [None, REQ]')
    p.add_option('--delim', default="\t", help="Delimiter used in the data file. [\\t]")
    p.add_option('--max', type="int", help="Size of output downsampled datasets. [None, REQ]")
    p.add_option('--min', type="int", default=0, help="Minimum dataset size for inclusion. Samples between min and max are included without downsampling [0]")
    p.add_option('--sd-mult', type="float", default=3.0, help="Standard deviation multiplier for zero-count filtering threshold [3.0]")
    p.add_option('--seed', type="int", default=None, help="Optional random seed for deterministic downsampling [None]")
    p.add_option('--failed-samples-output', default=None, help="Optional TSV output path for samples removed by this QC step.")
    p.add_option('--all-zero-counts-output', default=None, help="Optional TSV output path for zero counts from all considered samples.")
    p.add_option('--downsampled-counts-output', default=None, help="Optional TSV output path for full downsampled counts matrix before zero-count filtering.")
    p.add_option('--hist-output', default=None, help="Optional histogram output path (e.g., .png).")

    opts, args = p.parse_args()

    if not opts.data or not opts.out or opts.max is None:
        p.error("--data, --out, and --max are required.")
    if opts.seed is not None:
        random.seed(opts.seed)

    # Print command used
    print("Command run: '%s'" % ("  ".join(sys.argv)))
    
    # Read in data file
    if mt is not None:
        dataD = mt.parseCounts(opts.data, delim=opts.delim)
    else:
        print(
            "Warning: matrixtools not installed; using pandas fallback parser.",
            file=sys.stderr,
        )
        dataD = parse_counts_fallback(opts.data, delim=opts.delim)
    
    # Initiate dictionary to hold downsampled data
    toOut = {}
    
    # Step through each sample in the input data matrix to create
    # downsampled samples (or unchanged samples if within min/max).
    for s, dd in dataD.items():
        
        numReads = sum(dd.values())
        
        # If the total number of reads is greater than the max set by user, then downsample
        if numReads > opts.max:
            # Downsample without expanding reads into an in-memory list/dict.
            thisD = sampleCountsFromCounts(dd, opts.max)
            toOut[s] = thisD
    
        elif numReads >= opts.min:
            toOut[s] = dd

    if not toOut:
        raise ValueError("No samples met the --min threshold; nothing to output.")

    # Compute zero-counts on the downsampled matrix
    zeroCountD = {s: sum(1 for _, v in dd.items() if v == 0) for s, dd in toOut.items()}
    zeroCounts = list(zeroCountD.values())

    # Determine one-sided zero-count threshold
    zeroMean = statistics.mean(zeroCounts)
    # If only one sample survives min/max filtering, stdev is defined as 0.
    zeroStd = statistics.pstdev(zeroCounts) if len(zeroCounts) > 1 else 0.0
    upper = zeroMean + (opts.sd_mult * zeroStd)

    # Keep samples at or below the upper threshold.
    keepSamples = [s for s, z in zeroCountD.items() if z <= upper]
    if not keepSamples:
        raise ValueError("No samples pass the zero-count distribution filter.")
    keepSet = set(keepSamples)
    failedSamples = [s for s in toOut if s not in keepSet]

    # Output passing sample names only (no header), for downstream pepsirf subjoin.
    writePassSamples(opts.out, keepSamples)
    if opts.failed_samples_output:
        failedSampleZeroCounts = {s: zeroCountD[s] for s in failedSamples}
        writeFailedSamplesTsv(
            opts.failed_samples_output, failedSamples, STEP_NAME, failedSampleZeroCounts
        )
    if opts.all_zero_counts_output:
        writeAllZeroCountsTsv(opts.all_zero_counts_output, zeroCountD, keepSet)
    if opts.downsampled_counts_output:
        writeDownsampledCountsTsv(opts.downsampled_counts_output, toOut, delim=opts.delim)
    if opts.hist_output:
        writeZeroCountHistogram(opts.hist_output, zeroCounts, upper)

    # Summary
    print(f"Samples considered after min/max processing: {len(toOut)}")
    print(f"Zero-count mean: {zeroMean:.4f}")
    print(f"Zero-count std dev: {zeroStd:.4f}")
    print(f"Zero-count upper threshold: {upper:.4f}")
    print(f"Samples retained: {len(keepSamples)}")
    print(f"Samples failed {STEP_NAME}: {len(failedSamples)}")
    print(f"Pass-list output file: {opts.out}")
    if opts.hist_output:
        print(f"Histogram output: {opts.hist_output}")
    
def sampleCountsFromCounts(infoD: dict[str, Any], num2samp: int) -> dict[str, int]:
    """Downsample counts without materializing per-read expanded arrays."""
    population = []
    counts = []
    for each, count in infoD.items():
        c = int(count)
        if c > 0:
            population.append(each)
            counts.append(c)

    sampled = random.sample(population, k=num2samp, counts=counts)
    sampled_counts = Counter(sampled)

    out = {k: 0 for k in infoD}
    for k, v in sampled_counts.items():
        out[k] = v
    return out


def writeFailedSamplesTsv(
    path: str, failedSamples: list[str], stepName: str, sampleValues: dict[str, int]
) -> None:
    """Write per-sample QC failures and their failure values."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "Failed QC Step", "Failure Value"])
        for sample in failedSamples:
            writer.writerow([sample, stepName, sampleValues[sample]])


def writePassSamples(path: str, sampleNames: list[str]) -> None:
    """Write passing sample names, one per line, no header."""
    with open(path, "w", encoding="utf-8") as handle:
        for sample in sampleNames:
            handle.write(f"{sample}\n")


def writeAllZeroCountsTsv(
    path: str, sampleZeroCounts: dict[str, int], keepSet: set[str]
) -> None:
    """Write zero counts for all samples considered by this QC step."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "Zero Count", "Passed Filter"])
        for sample in sorted(sampleZeroCounts):
            writer.writerow([sample, sampleZeroCounts[sample], sample in keepSet])


def writeDownsampledCountsTsv(
    path: str, sampleCounts: dict[str, dict[str, int]], delim: str = "\t"
) -> None:
    """Write the full downsampled counts matrix (all considered samples)."""
    if not sampleCounts:
        raise ValueError("No downsampled counts available to write.")

    sampleNames = list(sampleCounts.keys())

    # Preserve the first sample's peptide order, then append any unseen peptides.
    firstSample = sampleNames[0]
    peptideOrder = list(sampleCounts[firstSample].keys())
    seen = set(peptideOrder)
    for sample in sampleNames[1:]:
        for peptide in sampleCounts[sample].keys():
            if peptide not in seen:
                peptideOrder.append(peptide)
                seen.add(peptide)

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter=delim)
        writer.writerow(["Sequence name", *sampleNames])
        for peptide in peptideOrder:
            writer.writerow([peptide, *[sampleCounts[s].get(peptide, 0) for s in sampleNames]])


def writeZeroCountHistogram(path: str, zeroCounts: list[int], upper: float) -> None:
    """Write a histogram of zero-count distribution and acceptance bounds."""
    plt.figure(figsize=(8, 5))
    plt.hist(zeroCounts, bins=30, edgecolor="black", alpha=0.8)
    plt.axvline(upper, color="red", linestyle="--", linewidth=1, label="Upper threshold")
    plt.xlabel("Zero counts per sample (downsampled)")
    plt.ylabel("Frequency")
    plt.title("Distribution of zero counts")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

if __name__ == "__main__":
    main()
