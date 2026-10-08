#!/usr/bin/env python3
"""Pass a count matrix through a disabled QC filter without computing metrics."""

import argparse
import csv
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--pass-output", required=True, type=Path)
    parser.add_argument("--failed-output", required=True, type=Path)
    parser.add_argument("--metrics-output", required=True, type=Path)
    parser.add_argument("--metric", required=True, choices=["read-counts", "zero-counts"])
    args = parser.parse_args()

    with args.input.open(newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    if not header:
        raise ValueError("Input count matrix has no header")
    samples = header[1:]
    shutil.copyfile(args.input, args.output)
    args.pass_output.write_text("".join(f"{sample}\n" for sample in samples))
    with args.failed_output.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerow(
            ["Sample", "Failed QC Step", "Failure Value"]
        )
    with args.metrics_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        if args.metric == "read-counts":
            writer.writerow(["Sample", "Total Read Counts"])
            writer.writerows([sample, ""] for sample in samples)
        else:
            writer.writerow(["Sample", "Zero Count", "Passed Filter"])
            writer.writerows([sample, "", ""] for sample in samples)
    print(f"{args.metric} QC disabled; retained all {len(samples)} samples. Metrics not computed.")


if __name__ == "__main__":
    main()
