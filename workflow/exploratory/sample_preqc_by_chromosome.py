#!/usr/bin/env python3

import argparse
import csv
import math
import statistics
from pathlib import Path

ROBUST_SCALE = 1.4826
PRIMARY_METRICS = ("average_depth", "missing_rate", "het_rate")
SECONDARY_METRICS = ("n_variants",)
ALL_Z_METRICS = PRIMARY_METRICS + SECONDARY_METRICS


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Summarize per-chromosome per-sample bcftools stats for diagnostic "
            "individual pre-QC. Flags are relative to the other samples within "
            "the same chromosome and never imply automatic exclusion."
        )
    )
    p.add_argument("--stats-dir", required=True)
    p.add_argument("--metrics-out", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--robust-z-threshold", type=float, default=3.5)
    return p.parse_args()


def safe_ratio(num, den):
    return float("nan") if den == 0 else num / den


def parse_stats(path, chrom):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("PSC\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 14:
                raise ValueError(f"Unexpected PSC record in {path}: {line.rstrip()}")

            sample = f[2]
            n_ref_hom = int(f[3])
            n_nonref_hom = int(f[4])
            n_hets = int(f[5])
            n_ts = int(f[6])
            n_tv = int(f[7])
            avg_dp = float(f[9])
            n_singletons = int(f[10])
            n_hap_ref = int(f[11])
            n_hap_alt = int(f[12])
            n_missing = int(f[13])

            n_called = n_ref_hom + n_nonref_hom + n_hets + n_hap_ref + n_hap_alt
            n_sites = n_called + n_missing
            n_variants = n_nonref_hom + n_hets + n_hap_alt

            rows.append({
                "chromosome": chrom,
                "sample": sample,
                "average_depth": avg_dp,
                "n_missing": n_missing,
                "missing_rate": safe_ratio(n_missing, n_sites),
                "n_hets": n_hets,
                "het_rate": safe_ratio(n_hets, n_called),
                "n_variants": n_variants,
                "n_singletons": n_singletons,
                "tstv": safe_ratio(n_ts, n_tv),
            })
    if not rows:
        raise ValueError(f"No PSC records found in {path}")
    return rows


def robust_center_scale(rows, metric):
    vals = [float(r[metric]) for r in rows if math.isfinite(float(r[metric]))]
    med = statistics.median(vals)
    mad = statistics.median(abs(x - med) for x in vals)
    sigma = ROBUST_SCALE * mad
    return med, sigma


def robust_z(value, med, sigma):
    if not math.isfinite(value):
        return float("nan")
    if sigma == 0:
        if value == med:
            return 0.0
        return math.copysign(float("inf"), value - med)
    return (value - med) / sigma


def fmt(x):
    if isinstance(x, float):
        if math.isnan(x):
            return "NA"
        if math.isinf(x):
            return "inf" if x > 0 else "-inf"
        return f"{x:.8g}"
    return str(x)


def main():
    args = parse_args()
    all_rows = []

    for chrom in range(1, 23):
        path = Path(args.stats_dir) / f"chr{chrom}.bcftools.stats.txt"
        chr_rows = parse_stats(path, chrom)

        stats = {m: robust_center_scale(chr_rows, m) for m in ALL_Z_METRICS}
        for row in chr_rows:
            primary_flags = []
            secondary_flags = []
            for metric in ALL_Z_METRICS:
                med, sigma = stats[metric]
                z = robust_z(float(row[metric]), med, sigma)
                row[f"robust_z_{metric}"] = z
                if math.isfinite(z) and abs(z) > args.robust_z_threshold:
                    if metric in PRIMARY_METRICS:
                        primary_flags.append(metric)
                    else:
                        secondary_flags.append(metric)
                elif not math.isfinite(z):
                    if metric in PRIMARY_METRICS:
                        primary_flags.append(metric)
                    else:
                        secondary_flags.append(metric)

            row["primary_flags"] = ",".join(primary_flags)
            row["secondary_flags"] = ",".join(secondary_flags)
            all_rows.append(row)

    metric_fields = [
        "chromosome", "sample",
        "average_depth", "missing_rate", "het_rate", "n_variants",
        "n_missing", "n_hets", "n_singletons", "tstv",
        "robust_z_average_depth", "robust_z_missing_rate",
        "robust_z_het_rate", "robust_z_n_variants",
        "primary_flags", "secondary_flags",
    ]

    Path(args.metrics_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.metrics_out, "w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=metric_fields, delimiter="\t")
        w.writeheader()
        for row in all_rows:
            w.writerow({k: fmt(row[k]) for k in metric_fields})

    samples = sorted({r["sample"] for r in all_rows})
    summary_fields = [
        "sample",
        "n_chromosomes",
        "chromosomes_flagged_any_primary",
        "chromosomes_flagged_average_depth",
        "chromosomes_flagged_missing_rate",
        "chromosomes_flagged_het_rate",
        "chromosomes_flagged_n_variants",
        "median_robust_z_average_depth",
        "median_robust_z_missing_rate",
        "median_robust_z_het_rate",
        "median_robust_z_n_variants",
    ]

    with open(args.summary_out, "w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        w.writeheader()

        for sample in samples:
            rows = [r for r in all_rows if r["sample"] == sample]

            def count_metric(metric):
                return sum(
                    metric in (r["primary_flags"] + "," + r["secondary_flags"]).split(",")
                    for r in rows
                )

            def median_z(metric):
                vals = [
                    float(r[f"robust_z_{metric}"])
                    for r in rows
                    if math.isfinite(float(r[f"robust_z_{metric}"]))
                ]
                return statistics.median(vals) if vals else float("nan")

            out = {
                "sample": sample,
                "n_chromosomes": len(rows),
                "chromosomes_flagged_any_primary": sum(bool(r["primary_flags"]) for r in rows),
                "chromosomes_flagged_average_depth": count_metric("average_depth"),
                "chromosomes_flagged_missing_rate": count_metric("missing_rate"),
                "chromosomes_flagged_het_rate": count_metric("het_rate"),
                "chromosomes_flagged_n_variants": count_metric("n_variants"),
                "median_robust_z_average_depth": median_z("average_depth"),
                "median_robust_z_missing_rate": median_z("missing_rate"),
                "median_robust_z_het_rate": median_z("het_rate"),
                "median_robust_z_n_variants": median_z("n_variants"),
            }
            w.writerow({k: fmt(out[k]) for k in summary_fields})


if __name__ == "__main__":
    main()
