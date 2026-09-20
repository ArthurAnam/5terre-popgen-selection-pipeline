#!/usr/bin/env python3

import argparse
import collections
import csv
import math
import random
from pathlib import Path

import matplotlib.pyplot as plt

NUMERIC_FIELDS = [
    "QUAL",
    "VQSLOD",
    "QD",
    "FS",
    "SOR",
    "MQ",
    "BaseQRankSum",
    "MQRankSum",
    "ReadPosRankSum",
    "INFO_DP",
    "ExcessHet",
    "InbreedingCoeff",
    "MAF",
    "F_MISSING",
]

PLOT_FIELDS = [
    "QUAL",
    "VQSLOD",
    "QD",
    "FS",
    "SOR",
    "MQ",
    "BaseQRankSum",
    "MQRankSum",
    "ReadPosRankSum",
    "INFO_DP",
    "ExcessHet",
    "InbreedingCoeff",
]


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Explore post-VQSR site annotations without applying additional hard filters. "
            "Exact counts are computed over all records; quantiles and plots use a fixed-seed "
            "reservoir sample."
        )
    )
    p.add_argument("--summary-out", required=True)
    p.add_argument("--inventory-out", required=True)
    p.add_argument("--culprit-out", required=True)
    p.add_argument("--plot-out", required=True)
    p.add_argument("--reservoir-size", type=int, default=100000)
    p.add_argument("--seed", type=int, default=20260918)
    return p.parse_args()


def parse_num(text):
    if text in ("", ".", "NA", None):
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


class Accumulator:
    def __init__(self, reservoir_size, rng):
        self.reservoir_size = reservoir_size
        self.rng = rng
        self.n_total = 0
        self.n_nonmissing = 0
        self.sum = 0.0
        self.min = None
        self.max = None
        self.sample = []

    def add(self, value):
        self.n_total += 1
        if value is None:
            return
        self.n_nonmissing += 1
        self.sum += value
        self.min = value if self.min is None else min(self.min, value)
        self.max = value if self.max is None else max(self.max, value)

        if len(self.sample) < self.reservoir_size:
            self.sample.append(value)
        else:
            j = self.rng.randrange(self.n_nonmissing)
            if j < self.reservoir_size:
                self.sample[j] = value

    @property
    def n_missing(self):
        return self.n_total - self.n_nonmissing

    @property
    def mean(self):
        return None if self.n_nonmissing == 0 else self.sum / self.n_nonmissing


def quantile(sorted_vals, p):
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = p * (len(sorted_vals) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def fmt(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.8g}"
    return str(value)


def is_snp(ref, alt):
    return len(ref) == 1 and len(alt) == 1 and "," not in alt


def is_palindromic(ref, alt):
    pair = {ref.upper(), alt.upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}


def main():
    args = parse_args()
    rng = random.Random(args.seed)
    acc = {field: Accumulator(args.reservoir_size, rng) for field in NUMERIC_FIELDS}

    counts = collections.Counter()
    culprit_counts = collections.Counter()

    reader = csv.DictReader(
        (line for line in __import__("sys").stdin if line.strip()),
        delimiter="\t",
    )

    required = {
        "CHROM", "POS", "REF", "ALT", "AC", "AN", "culprit",
        "POSITIVE_TRAIN_SITE", "NEGATIVE_TRAIN_SITE",
        *NUMERIC_FIELDS,
    }
    if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
        missing = sorted(required.difference(set(reader.fieldnames or [])))
        raise ValueError("Missing required input columns: " + ",".join(missing))

    for row in reader:
        counts["total_sites"] += 1

        ref = row["REF"]
        alt = row["ALT"]
        if is_snp(ref, alt):
            counts["biallelic_snps"] += 1
        else:
            if "," in alt:
                counts["multiallelic_sites"] += 1
            if len(ref) != 1 or any(len(a) != 1 for a in alt.split(",")):
                counts["non_snp_sites"] += 1

        if "," not in alt and is_palindromic(ref, alt):
            counts["palindromic_AT_CG"] += 1

        maf = parse_num(row["MAF"])
        if maf is not None and maf == 0:
            counts["monomorphic_after_sample_qc"] += 1

        f_missing = parse_num(row["F_MISSING"])
        if f_missing is not None and f_missing > 0.05:
            counts["site_missingness_gt_0.05"] += 1

        if row["POSITIVE_TRAIN_SITE"] not in ("", ".", "0"):
            counts["positive_train_sites"] += 1
        if row["NEGATIVE_TRAIN_SITE"] not in ("", ".", "0"):
            counts["negative_train_sites"] += 1

        culprit = row["culprit"] if row["culprit"] not in ("", ".") else "MISSING"
        culprit_counts[culprit] += 1

        for field in NUMERIC_FIELDS:
            acc[field].add(parse_num(row[field]))

    Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)

    quantiles = [
        ("q01", 0.01),
        ("q05", 0.05),
        ("q25", 0.25),
        ("median", 0.50),
        ("q75", 0.75),
        ("q95", 0.95),
        ("q99", 0.99),
    ]

    with open(args.summary_out, "w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "metric", "n_total", "n_nonmissing", "n_missing", "min",
            *[name for name, _ in quantiles],
            "max", "mean", "quantile_basis",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for metric in NUMERIC_FIELDS:
            a = acc[metric]
            sample_sorted = sorted(a.sample)
            out = {
                "metric": metric,
                "n_total": a.n_total,
                "n_nonmissing": a.n_nonmissing,
                "n_missing": a.n_missing,
                "min": fmt(a.min),
                "max": fmt(a.max),
                "mean": fmt(a.mean),
                "quantile_basis": f"fixed_seed_reservoir_n={len(a.sample)}",
            }
            for name, p in quantiles:
                out[name] = fmt(quantile(sample_sorted, p))
            writer.writerow(out)

    inventory_order = [
        "total_sites",
        "biallelic_snps",
        "multiallelic_sites",
        "non_snp_sites",
        "monomorphic_after_sample_qc",
        "palindromic_AT_CG",
        "site_missingness_gt_0.05",
        "positive_train_sites",
        "negative_train_sites",
    ]
    with open(args.inventory_out, "w", encoding="utf-8") as handle:
        handle.write("metric\tcount\n")
        for key in inventory_order:
            handle.write(f"{key}\t{counts[key]}\n")

    with open(args.culprit_out, "w", encoding="utf-8") as handle:
        handle.write("culprit\tcount\n")
        for key, value in culprit_counts.most_common():
            handle.write(f"{key}\t{value}\n")

    fig, axes = plt.subplots(3, 4, figsize=(15, 10))
    for ax, metric in zip(axes.flat, PLOT_FIELDS):
        values = sorted(acc[metric].sample)
        if not values:
            ax.set_title(metric + " (no data)")
            ax.axis("off")
            continue

        lo = quantile(values, 0.01)
        hi = quantile(values, 0.99)
        central = [v for v in values if lo <= v <= hi]
        if not central:
            central = values

        ax.hist(central, bins=60)
        ax.set_title(metric)
        ax.set_xlabel("value (1st-99th percentile shown)")
        ax.set_ylabel("sampled sites")

    fig.suptitle("Post-VQSR site annotation exploration (46 samples, autosomes 1-22)")
    fig.tight_layout()
    fig.savefig(args.plot_out, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
