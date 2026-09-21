#!/usr/bin/env python3

import argparse
import csv
import math
import statistics
from pathlib import Path


CORE_METRICS = ("n_singletons", "tstv", "average_depth", "n_variants")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Parse bcftools stats PSC records and reproduce the article-like "
            "sample QC based on five median absolute deviations."
        )
    )
    parser.add_argument("--stats", required=True)
    parser.add_argument("--metrics-out", required=True)
    parser.add_argument("--outliers-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--mad-threshold", type=float, default=5.0)
    return parser.parse_args()


def safe_ratio(num, den):
    return float("nan") if den == 0 else num / den


def parse_psc(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("PSC\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 14:
                raise ValueError(
                    f"Unexpected PSC record with {len(fields)} columns: {line.rstrip()}"
                )

            sample = fields[2]
            n_ref_hom = int(fields[3])
            n_nonref_hom = int(fields[4])
            n_hets = int(fields[5])
            n_transitions = int(fields[6])
            n_transversions = int(fields[7])
            n_indels = int(fields[8])
            average_depth = float(fields[9])
            n_singletons = int(fields[10])
            n_hap_ref = int(fields[11])
            n_hap_alt = int(fields[12])
            n_missing = int(fields[13])

            n_variants = n_nonref_hom + n_hets + n_hap_alt
            n_sites = (
                n_ref_hom
                + n_nonref_hom
                + n_hets
                + n_hap_ref
                + n_hap_alt
                + n_missing
            )

            rows.append(
                {
                    "sample": sample,
                    "n_ref_hom": n_ref_hom,
                    "n_nonref_hom": n_nonref_hom,
                    "n_hets": n_hets,
                    "n_transitions": n_transitions,
                    "n_transversions": n_transversions,
                    "tstv": safe_ratio(n_transitions, n_transversions),
                    "average_depth": average_depth,
                    "n_singletons": n_singletons,
                    "singleton_fraction": safe_ratio(n_singletons, n_variants),
                    "n_variants": n_variants,
                    "n_missing": n_missing,
                    "missing_rate": safe_ratio(n_missing, n_sites),
                    "n_indels": n_indels,
                    "n_hap_ref": n_hap_ref,
                    "n_hap_alt": n_hap_alt,
                }
            )

    if not rows:
        raise ValueError(
            "No PSC records found. Run bcftools stats with per-sample statistics enabled: -s -."
        )
    return rows


def finite_values(rows, key):
    vals = [float(row[key]) for row in rows]
    vals = [x for x in vals if math.isfinite(x)]
    if not vals:
        raise ValueError(f"No finite values available for metric {key}")
    return vals


def robust_summary(rows, key):
    vals = finite_values(rows, key)
    med = statistics.median(vals)
    mad = statistics.median(abs(x - med) for x in vals)
    return med, mad


def deviation(value, median, mad):
    if not math.isfinite(value):
        return float("inf")
    if mad == 0:
        return 0.0 if value == median else float("inf")
    return abs(value - median) / mad


def fmt(value):
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "inf"
        return f"{value:.8g}"
    return str(value)


def main():
    args = parse_args()
    rows = parse_psc(args.stats)

    metric_stats = {metric: robust_summary(rows, metric) for metric in CORE_METRICS}

    for row in rows:
        flagged = []
        for metric in CORE_METRICS:
            med, mad = metric_stats[metric]
            dev = deviation(float(row[metric]), med, mad)
            row[f"mad_dev_{metric}"] = dev
            if dev > args.mad_threshold:
                flagged.append(metric)
        row["outlier_metrics"] = ",".join(flagged)
        row["is_outlier"] = "YES" if flagged else "NO"

    fieldnames = [
        "sample",
        "n_variants",
        "n_singletons",
        "singleton_fraction",
        "tstv",
        "average_depth",
        "n_missing",
        "missing_rate",
        "n_ref_hom",
        "n_nonref_hom",
        "n_hets",
        "n_transitions",
        "n_transversions",
        "n_hap_ref",
        "n_hap_alt",
        "n_indels",
        "mad_dev_n_singletons",
        "mad_dev_tstv",
        "mad_dev_average_depth",
        "mad_dev_n_variants",
        "outlier_metrics",
        "is_outlier",
    ]

    Path(args.metrics_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.metrics_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row[k]) for k in fieldnames})

    outliers = [row for row in rows if row["is_outlier"] == "YES"]
    with open(args.outliers_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in outliers:
            writer.writerow({k: fmt(row[k]) for k in fieldnames})

    with open(args.summary_out, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"n_samples\t{len(rows)}\n")
        handle.write(f"mad_threshold\t{args.mad_threshold:g}\n")
        handle.write(f"n_outliers\t{len(outliers)}\n")
        handle.write(
            "outlier_samples\t"
            + (",".join(row["sample"] for row in outliers) if outliers else "NONE")
            + "\n"
        )
        for metric in CORE_METRICS:
            med, mad = metric_stats[metric]
            handle.write(f"median_{metric}\t{fmt(med)}\n")
            handle.write(f"mad_{metric}\t{fmt(mad)}\n")


if __name__ == "__main__":
    main()
