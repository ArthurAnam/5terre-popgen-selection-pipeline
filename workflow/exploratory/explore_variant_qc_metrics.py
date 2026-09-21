#!/usr/bin/env python3

import argparse
import csv
import math
import random
from pathlib import Path

import matplotlib.pyplot as plt

NUMERIC_FIELDS = ["QUAL", "MAF", "F_MISSING"]


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Describe current variant-QC metrics without applying filters. "
            "VQSR bookkeeping fields are deliberately excluded because VQSR "
            "was completed upstream."
        )
    )
    p.add_argument("--summary-out", required=True)
    p.add_argument("--inventory-out", required=True)
    p.add_argument("--plot-out", required=True)
    p.add_argument("--reservoir-size", type=int, default=100000)
    p.add_argument("--seed", type=int, default=20260920)
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


def is_biallelic_snp(ref, alt):
    return "," not in alt and len(ref) == 1 and len(alt) == 1


def is_palindromic(ref, alt):
    pair = {ref.upper(), alt.upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}


def log_edges(lo, hi, n_bins):
    lo = max(lo, 1e-9)
    if hi <= lo:
        hi = lo * 1.01
    a = math.log10(lo)
    b = math.log10(hi)
    return [10 ** (a + (b - a) * i / n_bins) for i in range(n_bins + 1)]


def main():
    args = parse_args()
    rng = random.Random(args.seed)
    acc = {field: Accumulator(args.reservoir_size, rng) for field in NUMERIC_FIELDS}

    counts = {
        "total_sites": 0,
        "biallelic_snps": 0,
        "multiallelic_sites": 0,
        "non_snp_sites": 0,
        "monomorphic_after_sample_qc": 0,
        "palindromic_AT_CG": 0,
        "site_missingness_gt_0.05": 0,
    }

    reader = csv.DictReader(
        (line for line in __import__("sys").stdin if line.strip()),
        delimiter="\t",
    )
    required = {"CHROM", "POS", "REF", "ALT", "QUAL", "MAF", "F_MISSING"}
    if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
        missing = sorted(required.difference(set(reader.fieldnames or [])))
        raise ValueError("Missing required input columns: " + ",".join(missing))

    for row in reader:
        counts["total_sites"] += 1
        ref = row["REF"]
        alt = row["ALT"]

        if is_biallelic_snp(ref, alt):
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

        for field in NUMERIC_FIELDS:
            acc[field].add(parse_num(row[field]))

    Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)

    quantiles = [
        ("q01", 0.01), ("q05", 0.05), ("q25", 0.25), ("median", 0.50),
        ("q75", 0.75), ("q95", 0.95), ("q99", 0.99),
    ]

    with open(args.summary_out, "w", newline="", encoding="utf-8") as handle:
        fields = [
            "metric", "n_total", "n_nonmissing", "n_missing", "min",
            *[name for name, _ in quantiles], "max", "mean", "quantile_basis",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for metric in NUMERIC_FIELDS:
            a = acc[metric]
            vals = sorted(a.sample)
            row = {
                "metric": metric,
                "n_total": a.n_total,
                "n_nonmissing": a.n_nonmissing,
                "n_missing": a.n_missing,
                "min": fmt(a.min),
                "max": fmt(a.max),
                "mean": fmt(a.mean),
                "quantile_basis": f"fixed_seed_reservoir_n={len(vals)}",
            }
            for name, p in quantiles:
                row[name] = fmt(quantile(vals, p))
            writer.writerow(row)

    with open(args.inventory_out, "w", encoding="utf-8") as handle:
        handle.write("metric\tcount\n")
        for key, value in counts.items():
            handle.write(f"{key}\t{value}\n")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # QUAL: log-scaled x axis so the long right tail does not hide the
    # informative low-to-middle range. Tick values remain QUAL values.
    qual = sorted(acc["QUAL"].sample)
    q01 = quantile(qual, 0.01)
    q50 = quantile(qual, 0.50)
    q99 = quantile(qual, 0.99)
    qual_central = [v for v in qual if q01 <= v <= q99 and v > 0]
    axes[0].hist(qual_central, bins=log_edges(q01, q99, 60))
    axes[0].set_xscale("log")
    axes[0].axvline(q50, linestyle="--")
    axes[0].set_title("QUAL (1st-99th percentile)")
    axes[0].set_xlabel("QUAL, log-scaled axis")
    axes[0].set_ylabel("sampled sites")
    axes[0].text(
        0.98, 0.96,
        f"min={fmt(acc['QUAL'].min)}\nmedian={fmt(q50)}\nq99={fmt(q99)}\nmax={fmt(acc['QUAL'].max)}",
        transform=axes[0].transAxes,
        ha="right", va="top",
    )

    # MAF: QC-relevant zoom on the rare/low-frequency end. MAF is descriptive
    # here; no global MAF hard filter is part of core QC.
    maf = acc["MAF"].sample
    maf_zoom = [v for v in maf if 0 <= v <= 0.10]
    axes[1].hist(maf_zoom, bins=40, range=(0, 0.10))
    axes[1].set_xlim(0, 0.10)
    axes[1].set_title("MAF (zoom: 0-0.10)")
    axes[1].set_xlabel("minor allele frequency")
    axes[1].set_ylabel("sampled sites")
    axes[1].text(
        0.98, 0.96,
        "descriptive only\nno MAF hard filter",
        transform=axes[1].transAxes,
        ha="right", va="top",
    )

    # Missingness: show the range around the planned 5% threshold rather than
    # truncating at q99, which hid the QC-relevant tail in the previous plot.
    miss = acc["F_MISSING"].sample
    miss_zoom = [v for v in miss if 0 <= v <= 0.10]
    axes[2].hist(miss_zoom, bins=46, range=(0, 0.10))
    axes[2].set_yscale("log")
    axes[2].axvline(0.05, linestyle="--")
    axes[2].set_xlim(0, 0.10)
    axes[2].set_title("Variant missingness (0-10%)")
    axes[2].set_xlabel("fraction of missing genotypes")
    axes[2].set_ylabel("sampled sites, log scale")
    axes[2].text(
        0.98, 0.96,
        f"planned cutoff: >5%\nexact sites >5%: {counts['site_missingness_gt_0.05']}",
        transform=axes[2].transAxes,
        ha="right", va="top",
    )

    fig.suptitle("Variant-QC exploration: descriptive only, no filter applied")
    fig.tight_layout()
    fig.savefig(args.plot_out, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
