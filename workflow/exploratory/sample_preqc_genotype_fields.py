#!/usr/bin/env python3

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

ROBUST_SCALE = 1.4826

Z_METRICS = (
    "mean_dp",
    "mean_gq",
    "mean_het_alt_fraction",
    "mean_het_abs_balance_deviation",
    "mean_hom_ref_alt_fraction",
    "mean_hom_alt_ref_fraction",
)


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Summarize GT/DP/GQ/AD over a deterministic genome-wide subset of "
            "autosomal sites for diagnostic per-sample pre-QC. No hard DP, GQ "
            "or allele-balance thresholds are applied."
        )
    )
    p.add_argument("--samples", required=True)
    p.add_argument("--metrics-out", required=True)
    return p.parse_args()


def safe_float(text):
    if text in ("", ".", None):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_ad(text):
    if text in ("", ".", None):
        return None
    parts = text.split(",")
    if len(parts) < 2:
        return None
    try:
        ref = int(parts[0])
        alt = int(parts[1])
    except ValueError:
        return None
    total = ref + alt
    if total <= 0:
        return None
    return ref, alt, total


def genotype_class(gt):
    if not gt or gt in (".", "./.", ".|."):
        return "missing"
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or "." in alleles:
        return "other"
    if alleles == ["0", "0"]:
        return "hom_ref"
    if alleles == ["1", "1"]:
        return "hom_alt"
    if set(alleles) == {"0", "1"}:
        return "het"
    return "other"


def mean_or_nan(total, n):
    return float("nan") if n == 0 else total / n


def robust_z(values):
    finite = [x for x in values if math.isfinite(x)]
    if not finite:
        return [float("nan")] * len(values)
    med = statistics.median(finite)
    mad = statistics.median(abs(x - med) for x in finite)
    sigma = ROBUST_SCALE * mad

    out = []
    for x in values:
        if not math.isfinite(x):
            out.append(float("nan"))
        elif sigma == 0:
            out.append(0.0 if x == med else math.copysign(float("inf"), x - med))
        else:
            out.append((x - med) / sigma)
    return out


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

    with open(args.samples, "r", encoding="utf-8") as handle:
        samples = [line.strip() for line in handle if line.strip()]

    if not samples:
        raise ValueError("No samples found in sample list")

    rows = []
    for sample in samples:
        rows.append({
            "sample": sample,
            "sampled_sites": 0,
            "called_genotypes": 0,
            "het_genotypes": 0,
            "hom_ref_genotypes": 0,
            "hom_alt_genotypes": 0,
            "dp_available": 0,
            "sum_dp": 0.0,
            "gq_available": 0,
            "sum_gq": 0.0,
            "het_ad_available": 0,
            "sum_het_alt_fraction": 0.0,
            "sum_het_abs_balance_deviation": 0.0,
            "hom_ref_ad_available": 0,
            "sum_hom_ref_alt_fraction": 0.0,
            "hom_alt_ad_available": 0,
            "sum_hom_alt_ref_fraction": 0.0,
        })

    expected_fields = 2 + 4 * len(samples)
    n_lines = 0

    for line in sys.stdin:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != expected_fields:
            raise ValueError(
                f"Expected {expected_fields} fields per site, found {len(fields)} "
                f"at input line {n_lines + 1}"
            )
        n_lines += 1

        for i, row in enumerate(rows):
            base = 2 + 4 * i
            gt = fields[base]
            dp = safe_float(fields[base + 1])
            gq = safe_float(fields[base + 2])
            ad = parse_ad(fields[base + 3])

            row["sampled_sites"] += 1
            gt_class = genotype_class(gt)

            if gt_class != "missing":
                row["called_genotypes"] += 1

            if dp is not None:
                row["dp_available"] += 1
                row["sum_dp"] += dp

            if gq is not None:
                row["gq_available"] += 1
                row["sum_gq"] += gq

            if gt_class == "het":
                row["het_genotypes"] += 1
                if ad is not None:
                    ref, alt, total = ad
                    ab = alt / total
                    row["het_ad_available"] += 1
                    row["sum_het_alt_fraction"] += ab
                    row["sum_het_abs_balance_deviation"] += abs(ab - 0.5)

            elif gt_class == "hom_ref":
                row["hom_ref_genotypes"] += 1
                if ad is not None:
                    ref, alt, total = ad
                    row["hom_ref_ad_available"] += 1
                    row["sum_hom_ref_alt_fraction"] += alt / total

            elif gt_class == "hom_alt":
                row["hom_alt_genotypes"] += 1
                if ad is not None:
                    ref, alt, total = ad
                    row["hom_alt_ad_available"] += 1
                    row["sum_hom_alt_ref_fraction"] += ref / total

    if n_lines == 0:
        raise ValueError("No genotype-field records were received on stdin")

    for row in rows:
        row["call_rate_sampled"] = mean_or_nan(row["called_genotypes"], row["sampled_sites"])
        row["mean_dp"] = mean_or_nan(row["sum_dp"], row["dp_available"])
        row["mean_gq"] = mean_or_nan(row["sum_gq"], row["gq_available"])
        row["mean_het_alt_fraction"] = mean_or_nan(
            row["sum_het_alt_fraction"], row["het_ad_available"]
        )
        row["mean_het_abs_balance_deviation"] = mean_or_nan(
            row["sum_het_abs_balance_deviation"], row["het_ad_available"]
        )
        row["mean_hom_ref_alt_fraction"] = mean_or_nan(
            row["sum_hom_ref_alt_fraction"], row["hom_ref_ad_available"]
        )
        row["mean_hom_alt_ref_fraction"] = mean_or_nan(
            row["sum_hom_alt_ref_fraction"], row["hom_alt_ad_available"]
        )

    for metric in Z_METRICS:
        zvals = robust_z([float(row[metric]) for row in rows])
        for row, z in zip(rows, zvals):
            row[f"robust_z_{metric}"] = z

    fieldnames = [
        "sample",
        "sampled_sites",
        "called_genotypes",
        "call_rate_sampled",
        "het_genotypes",
        "hom_ref_genotypes",
        "hom_alt_genotypes",
        "dp_available",
        "mean_dp",
        "gq_available",
        "mean_gq",
        "het_ad_available",
        "mean_het_alt_fraction",
        "mean_het_abs_balance_deviation",
        "hom_ref_ad_available",
        "mean_hom_ref_alt_fraction",
        "hom_alt_ad_available",
        "mean_hom_alt_ref_fraction",
        "robust_z_mean_dp",
        "robust_z_mean_gq",
        "robust_z_mean_het_alt_fraction",
        "robust_z_mean_het_abs_balance_deviation",
        "robust_z_mean_hom_ref_alt_fraction",
        "robust_z_mean_hom_alt_ref_fraction",
    ]

    Path(args.metrics_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.metrics_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row[key]) for key in fieldnames})


if __name__ == "__main__":
    main()
