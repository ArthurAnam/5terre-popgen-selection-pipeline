#!/usr/bin/env python3

import argparse
import csv
import math
import statistics
from pathlib import Path

ROBUST_SCALE = 1.4826

PRIMARY_METRICS = {
    "average_depth": "two_sided",
    "missing_rate": "upper",
    "het_rate": "two_sided",
}

SECONDARY_METRICS = {
    "n_variants": "two_sided",
    "singleton_fraction": "two_sided",
    "tstv": "two_sided",
}

ALL_FLAG_METRICS = {**PRIMARY_METRICS, **SECONDARY_METRICS}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize autosomal per-sample bcftools stats for diagnostic sample pre-QC. "
            "Robust-z flags are for review only and never imply automatic exclusion."
        )
    )
    parser.add_argument("--stats", required=True)
    parser.add_argument("--metrics-out", required=True)
    parser.add_argument("--flags-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--robust-z-threshold", type=float, default=3.5)
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

            # This stage is restricted to autosomes, so diploid genotype counts
            # are the quantities used for call rate and heterozygosity.
            n_called = n_ref_hom + n_nonref_hom + n_hets
            n_sites = n_called + n_missing
            n_variants = n_nonref_hom + n_hets

            rows.append(
                {
                    "sample": sample,
                    "n_sites": n_sites,
                    "n_called": n_called,
                    "n_missing": n_missing,
                    "missing_rate": safe_ratio(n_missing, n_sites),
                    "call_rate": safe_ratio(n_called, n_sites),
                    "average_depth": average_depth,
                    "n_ref_hom": n_ref_hom,
                    "n_nonref_hom": n_nonref_hom,
                    "n_hets": n_hets,
                    "het_rate": safe_ratio(n_hets, n_called),
                    "n_variants": n_variants,
                    "n_singletons": n_singletons,
                    "singleton_fraction": safe_ratio(n_singletons, n_variants),
                    "n_transitions": n_transitions,
                    "n_transversions": n_transversions,
                    "tstv": safe_ratio(n_transitions, n_transversions),
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
    robust_sigma = ROBUST_SCALE * mad
    return med, mad, robust_sigma


def robust_z(value, median, robust_sigma):
    if not math.isfinite(value):
        return float("nan")
    if robust_sigma == 0:
        return 0.0 if value == median else math.copysign(float("inf"), value - median)
    return (value - median) / robust_sigma


def should_flag(z, direction, threshold):
    if not math.isfinite(z):
        return True
    if direction == "upper":
        return z > threshold
    if direction == "lower":
        return z < -threshold
    return abs(z) > threshold


def fmt(value):
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        return f"{value:.8g}"
    return str(value)


def main():
    args = parse_args()
    rows = parse_psc(args.stats)

    metric_stats = {
        metric: robust_summary(rows, metric) for metric in ALL_FLAG_METRICS
    }

    for row in rows:
        primary_flags = []
        secondary_flags = []

        for metric, direction in ALL_FLAG_METRICS.items():
            med, mad, robust_sigma = metric_stats[metric]
            z = robust_z(float(row[metric]), med, robust_sigma)
            row[f"robust_z_{metric}"] = z

            if should_flag(z, direction, args.robust_z_threshold):
                if metric in PRIMARY_METRICS:
                    primary_flags.append(metric)
                else:
                    secondary_flags.append(metric)

        row["primary_flags"] = ",".join(primary_flags)
        row["secondary_flags"] = ",".join(secondary_flags)
        row["n_primary_flags"] = len(primary_flags)
        row["n_secondary_flags"] = len(secondary_flags)
        row["review_flag"] = "YES" if primary_flags or secondary_flags else "NO"

    fieldnames = [
        "sample",
        "n_sites",
        "n_called",
        "n_missing",
        "missing_rate",
        "call_rate",
        "average_depth",
        "n_ref_hom",
        "n_nonref_hom",
        "n_hets",
        "het_rate",
        "n_variants",
        "n_singletons",
        "singleton_fraction",
        "n_transitions",
        "n_transversions",
        "tstv",
        "n_indels",
        "robust_z_average_depth",
        "robust_z_missing_rate",
        "robust_z_het_rate",
        "robust_z_n_variants",
        "robust_z_singleton_fraction",
        "robust_z_tstv",
        "primary_flags",
        "secondary_flags",
        "n_primary_flags",
        "n_secondary_flags",
        "review_flag",
    ]

    Path(args.metrics_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.metrics_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row[key]) for key in fieldnames})

    flagged_rows = [row for row in rows if row["review_flag"] == "YES"]
    with open(args.flags_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in flagged_rows:
            writer.writerow({key: fmt(row[key]) for key in fieldnames})

    primary_flagged = [row for row in rows if row["n_primary_flags"] > 0]
    secondary_only = [
        row
        for row in rows
        if row["n_primary_flags"] == 0 and row["n_secondary_flags"] > 0
    ]

    with open(args.summary_out, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"scope\tautosomes_1_22\n")
        handle.write(f"n_samples\t{len(rows)}\n")
        handle.write(f"robust_z_flag_threshold\t{args.robust_z_threshold:g}\n")
        handle.write("automatic_sample_exclusion\tNO\n")
        handle.write(f"n_review_flagged_samples\t{len(flagged_rows)}\n")
        handle.write(f"n_primary_flagged_samples\t{len(primary_flagged)}\n")
        handle.write(f"n_secondary_only_flagged_samples\t{len(secondary_only)}\n")
        handle.write(
            "review_flagged_samples\t"
            + (",".join(row["sample"] for row in flagged_rows) if flagged_rows else "NONE")
            + "\n"
        )
        handle.write(
            "primary_flagged_samples\t"
            + (",".join(row["sample"] for row in primary_flagged) if primary_flagged else "NONE")
            + "\n"
        )
        handle.write(
            "secondary_only_flagged_samples\t"
            + (",".join(row["sample"] for row in secondary_only) if secondary_only else "NONE")
            + "\n"
        )

        for metric in ALL_FLAG_METRICS:
            med, mad, robust_sigma = metric_stats[metric]
            handle.write(f"median_{metric}\t{fmt(med)}\n")
            handle.write(f"mad_{metric}\t{fmt(mad)}\n")
            handle.write(f"robust_sigma_{metric}\t{fmt(robust_sigma)}\n")


if __name__ == "__main__":
    main()
