#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Evaluate chromosome-specific Bonferroni HWE thresholds. "
            "Threshold per chromosome = alpha / number of tested variants."
        )
    )
    p.add_argument("--thresholds", required=True)
    p.add_argument("--failed-out", required=True)
    p.add_argument("--summary-out", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    threshold = {}
    expected = {}
    alpha = {}
    with open(args.thresholds, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = row["chromosome"]
            expected[chrom] = int(row["n_tested"])
            alpha[chrom] = float(row["alpha"])
            threshold[chrom] = float(row["bonferroni_threshold"])

    seen = defaultdict(int)
    failed = defaultdict(int)
    missing_hwe = defaultdict(int)

    reader = csv.DictReader(sys.stdin, delimiter="\t")
    required = {"CHROM", "POS", "REF", "ALT", "HWE"}
    if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
        raise ValueError("Missing required HWE input fields")

    with open(args.failed_out, "w", encoding="utf-8") as failed_handle:
        failed_handle.write("chromosome\tposition\tref\talt\thwe_p\tbonferroni_threshold\n")

        for row in reader:
            chrom = row["CHROM"]
            if chrom not in threshold:
                raise ValueError(f"No Bonferroni threshold defined for chromosome {chrom}")

            seen[chrom] += 1
            hwe_text = row["HWE"]
            if hwe_text in ("", ".", "NA"):
                missing_hwe[chrom] += 1
                continue

            p = float(hwe_text)
            if p < threshold[chrom]:
                failed[chrom] += 1
                failed_handle.write(
                    f"{chrom}\t{row['POS']}\t{row['REF']}\t{row['ALT']}\t"
                    f"{p:.12g}\t{threshold[chrom]:.12g}\n"
                )

    for chrom, n_expected in expected.items():
        if seen[chrom] != n_expected:
            raise ValueError(
                f"Chromosome {chrom}: expected {n_expected} records, observed {seen[chrom]}"
            )
        if missing_hwe[chrom] != 0:
            raise ValueError(
                f"Chromosome {chrom}: HWE was missing for {missing_hwe[chrom]} records"
            )

    with open(args.summary_out, "w", encoding="utf-8") as handle:
        handle.write(
            "chromosome\tn_tested\talpha\tbonferroni_threshold\t"
            "n_removed_hwe\tn_retained_hwe\n"
        )
        for chrom in sorted(expected, key=lambda x: int(x)):
            n = expected[chrom]
            n_fail = failed[chrom]
            handle.write(
                f"{chrom}\t{n}\t{alpha[chrom]:.12g}\t{threshold[chrom]:.12g}\t"
                f"{n_fail}\t{n - n_fail}\n"
            )


if __name__ == "__main__":
    main()
