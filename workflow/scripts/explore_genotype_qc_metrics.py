#!/usr/bin/env python3

import argparse
import csv
import math
import random
import sys
from pathlib import Path

import matplotlib.pyplot as plt

METRICS = ["DP", "GQ", "HET_ALT_FRACTION", "HET_ABS_BALANCE_DEVIATION"]


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Explore FORMAT-level DP, GQ and AD on a deterministic genome-wide "
            "subset. This is descriptive only and does not mask or filter genotypes."
        )
    )
    p.add_argument("--summary-out", required=True)
    p.add_argument("--plot-out", required=True)
    p.add_argument("--reservoir-size", type=int, default=200000)
    p.add_argument("--seed", type=int, default=20260920)
    return p.parse_args()


def parse_num(text):
    if text in ("", ".", None):
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def gt_class(gt):
    if gt in ("", ".", "./.", ".|.", None) or "." in gt:
        return "missing"
    a = gt.replace("|", "/").split("/")
    if len(a) != 2:
        return "other"
    if set(a) == {"0", "1"}:
        return "het"
    return "called"


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


class Reservoir:
    def __init__(self, size, rng):
        self.size = size
        self.rng = rng
        self.n = 0
        self.sum = 0.0
        self.min = None
        self.max = None
        self.values = []

    def add(self, value):
        if value is None:
            return
        self.n += 1
        self.sum += value
        self.min = value if self.min is None else min(self.min, value)
        self.max = value if self.max is None else max(self.max, value)
        if len(self.values) < self.size:
            self.values.append(value)
        else:
            j = self.rng.randrange(self.n)
            if j < self.size:
                self.values[j] = value

    @property
    def mean(self):
        return None if self.n == 0 else self.sum / self.n


def quantile(vals, p):
    if not vals:
        return None
    vals = sorted(vals)
    if len(vals) == 1:
        return vals[0]
    pos = p * (len(vals) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    frac = pos - lo
    return vals[lo] * (1 - frac) + vals[hi] * frac


def fmt(x):
    if x is None:
        return "NA"
    return f"{x:.8g}" if isinstance(x, float) else str(x)


def main():
    args = parse_args()
    rng = random.Random(args.seed)
    acc = {m: Reservoir(args.reservoir_size, rng) for m in METRICS}

    header = sys.stdin.readline().rstrip("\n").split("\t")
    if not header or header[0] != "CHROM":
        raise ValueError("Expected a tab-delimited header starting with CHROM")

    n_genotype_fields = len(header) - 2
    if n_genotype_fields % 4 != 0:
        raise ValueError("Expected GT,DP,GQ,AD repeated for each sample")
    n_samples = n_genotype_fields // 4

    n_sites = 0
    n_genotypes = 0
    n_called = 0
    n_missing = 0

    for line in sys.stdin:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != len(header):
            raise ValueError("Unexpected number of columns in genotype query")
        n_sites += 1

        for i in range(n_samples):
            base = 2 + 4 * i
            gt = fields[base]
            dp = parse_num(fields[base + 1])
            gq = parse_num(fields[base + 2])
            ad = parse_ad(fields[base + 3])

            n_genotypes += 1
            cls = gt_class(gt)
            if cls == "missing":
                n_missing += 1
            else:
                n_called += 1

            acc["DP"].add(dp)
            acc["GQ"].add(gq)

            if cls == "het" and ad is not None:
                ref, alt, total = ad
                ab = alt / total
                acc["HET_ALT_FRACTION"].add(ab)
                acc["HET_ABS_BALANCE_DEVIATION"].add(abs(ab - 0.5))

    if n_sites == 0:
        raise ValueError("No records received")

    Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)
    qs = [("q01", 0.01), ("q05", 0.05), ("q25", 0.25), ("median", 0.50),
          ("q75", 0.75), ("q95", 0.95), ("q99", 0.99)]

    with open(args.summary_out, "w", newline="", encoding="utf-8") as handle:
        fields = ["metric", "n_observations", "min", *[x[0] for x in qs],
                  "max", "mean", "quantile_basis"]
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for metric in METRICS:
            a = acc[metric]
            row = {
                "metric": metric,
                "n_observations": a.n,
                "min": fmt(a.min),
                "max": fmt(a.max),
                "mean": fmt(a.mean),
                "quantile_basis": f"fixed_seed_reservoir_n={len(a.values)}",
            }
            for name, p in qs:
                row[name] = fmt(quantile(a.values, p))
            w.writerow(row)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    for ax, metric, title, xlabel in [
        (axes[0], "DP", "Genotype depth (DP)", "FORMAT/DP"),
        (axes[1], "GQ", "Genotype quality (GQ)", "FORMAT/GQ"),
    ]:
        vals = sorted(acc[metric].values)
        lo = quantile(vals, 0.01)
        hi = quantile(vals, 0.99)
        central = [v for v in vals if lo <= v <= hi]
        ax.hist(central, bins=60)
        ax.axvline(quantile(vals, 0.50), linestyle="--")
        ax.set_title(title)
        ax.set_xlabel(f"{xlabel} (1st-99th percentile)")
        ax.set_ylabel("sampled genotypes")

    ab = acc["HET_ALT_FRACTION"].values
    axes[2].hist(ab, bins=50, range=(0, 1))
    axes[2].axvline(0.5, linestyle="--")
    axes[2].set_xlim(0, 1)
    axes[2].set_title("Allelic balance in heterozygotes")
    axes[2].set_xlabel("ALT / (REF + ALT) from FORMAT/AD")
    axes[2].set_ylabel("sampled heterozygous genotypes")

    fig.suptitle(
        f"Genotype-QC exploration: {n_samples} samples, {n_sites} sampled sites; no filter applied"
    )
    fig.tight_layout()
    fig.savefig(args.plot_out, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
