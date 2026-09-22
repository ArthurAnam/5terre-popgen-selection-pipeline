#!/usr/bin/env python3
import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

POP_ORDER = ["CT", "CEU", "FIN", "GBR", "IBS", "TSI"]
REFERENCE_POPS = ["CEU", "FIN", "GBR", "IBS", "TSI"]
ENDPOINTS = [
    ("FROH_GE1_5MB", "primary_froh_ge1.5mb"),
    ("FROH_GE5MB", "secondary_froh_ge5mb"),
]

def mean(x):
    return statistics.mean(x) if x else float("nan")

def median(x):
    return statistics.median(x) if x else float("nan")

def sd(x):
    return statistics.stdev(x) if len(x) > 1 else 0.0

def quantile(values, q):
    vals = sorted(values)
    if not vals:
        return float("nan")
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    frac = pos - lo
    return vals[lo] * (1.0 - frac) + vals[hi] * frac

def average_ranks(values):
    indexed = sorted(enumerate(values), key=lambda z: z[1])
    ranks = [0.0] * len(values)
    tie_sizes = []
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[indexed[k][0]] = rank
        tie_sizes.append(j - i)
        i = j
    return ranks, tie_sizes

def mann_whitney(x, y):
    n1, n2 = len(x), len(y)
    pooled = x + y
    ranks, ties = average_ranks(pooled)
    r1 = sum(ranks[:n1])
    u1 = r1 - n1 * (n1 + 1) / 2.0
    mean_u = n1 * n2 / 2.0
    n = n1 + n2
    tie_term = sum(t**3 - t for t in ties)
    variance = n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1))) if n > 1 else 0.0
    if variance <= 0:
        p = 1.0
        z = 0.0
    else:
        z = (u1 - mean_u) / math.sqrt(variance)
        p = math.erfc(abs(z) / math.sqrt(2.0))
    cliffs_delta = 2.0 * u1 / (n1 * n2) - 1.0
    return u1, z, p, cliffs_delta

def holm_adjust(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adjusted = [0.0] * m
    running = 0.0
    for rank, idx in enumerate(order):
        value = (m - rank) * pvals[idx]
        running = max(running, value)
        adjusted[idx] = min(1.0, running)
    return adjusted

def read_rows(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"population","IID","N_ROH_GE1_5MB","TOTAL_ROH_MB_GE1_5MB","FROH_GE1_5MB","MAX_ROH_MB","FROH_GE5MB"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit("ERROR: production individual table is missing required columns")
        for row in reader:
            row["N_ROH_GE1_5MB"] = int(row["N_ROH_GE1_5MB"])
            for field in ["TOTAL_ROH_MB_GE1_5MB","FROH_GE1_5MB","MAX_ROH_MB","FROH_GE5MB"]:
                row[field] = float(row[field])
            rows.append(row)
    return rows

def jitter(i, j):
    # Deterministic small offset without importing random.
    return i + (((j * 37) % 101) / 100.0 - 0.5) * 0.28

def make_box_scatter(by_pop, field, ylabel, png, pdf):
    data = [[r[field] for r in by_pop[p]] for p in POP_ORDER]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.boxplot(data, tick_labels=POP_ORDER, showfliers=False)
    for i, pop in enumerate(POP_ORDER, start=1):
        vals = data[i - 1]
        ax.scatter([jitter(i, j) for j in range(len(vals))], vals, s=13, alpha=0.55)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Population")
    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

def make_scatter(by_pop, png, pdf):
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for pop in POP_ORDER:
        x = [r["N_ROH_GE1_5MB"] for r in by_pop[pop]]
        y = [r["TOTAL_ROH_MB_GE1_5MB"] for r in by_pop[pop]]
        ax.scatter(x, y, s=24, alpha=0.65, label=pop)
    ax.set_xlabel("Number of ROH >=1.5 Mb")
    ax.set_ylabel("Total ROH length >=1.5 Mb (Mb)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--individual", required=True)
    p.add_argument("--descriptive-out", required=True)
    p.add_argument("--tests-out", required=True)
    p.add_argument("--froh-png", required=True)
    p.add_argument("--froh-pdf", required=True)
    p.add_argument("--long-png", required=True)
    p.add_argument("--long-pdf", required=True)
    p.add_argument("--scatter-png", required=True)
    p.add_argument("--scatter-pdf", required=True)
    p.add_argument("--readme-out", required=True)
    a = p.parse_args()

    rows = read_rows(a.individual)
    by_pop = defaultdict(list)
    for row in rows:
        by_pop[row["population"]].append(row)
    missing = [p for p in POP_ORDER if p not in by_pop]
    if missing:
        raise SystemExit(f"ERROR: missing populations: {missing}")

    Path(a.descriptive_out).parent.mkdir(parents=True, exist_ok=True)

    desc_fields = ["population","endpoint","n","mean","sd","q1","median","q3","min","max","zero_fraction"]
    with open(a.descriptive_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=desc_fields)
        writer.writeheader()
        for pop in POP_ORDER:
            for field, endpoint in ENDPOINTS:
                vals = [r[field] for r in by_pop[pop]]
                writer.writerow({
                    "population": pop, "endpoint": endpoint, "n": len(vals),
                    "mean": mean(vals), "sd": sd(vals), "q1": quantile(vals, 0.25),
                    "median": median(vals), "q3": quantile(vals, 0.75),
                    "min": min(vals), "max": max(vals),
                    "zero_fraction": sum(v == 0 for v in vals) / len(vals),
                })

    test_rows = []
    for field, endpoint in ENDPOINTS:
        raw_ps = []
        endpoint_rows = []
        x = [r[field] for r in by_pop["CT"]]
        for ref in REFERENCE_POPS:
            y = [r[field] for r in by_pop[ref]]
            u, z, pval, delta = mann_whitney(x, y)
            rec = {
                "endpoint": endpoint, "comparison": f"CT_vs_{ref}",
                "n_CT": len(x), "n_reference": len(y),
                "mean_CT": mean(x), "mean_reference": mean(y),
                "median_CT": median(x), "median_reference": median(y),
                "mean_difference_CT_minus_reference": mean(x) - mean(y),
                "median_difference_CT_minus_reference": median(x) - median(y),
                "mann_whitney_U_CT": u, "z_tie_corrected": z,
                "p_two_sided": pval, "cliffs_delta_CT_minus_reference": delta,
            }
            raw_ps.append(pval)
            endpoint_rows.append(rec)
        adjusted = holm_adjust(raw_ps)
        for rec, padj in zip(endpoint_rows, adjusted):
            rec["p_holm_within_endpoint"] = padj
            test_rows.append(rec)

    test_fields = ["endpoint","comparison","n_CT","n_reference","mean_CT","mean_reference",
                   "median_CT","median_reference","mean_difference_CT_minus_reference",
                   "median_difference_CT_minus_reference","mann_whitney_U_CT","z_tie_corrected",
                   "p_two_sided","p_holm_within_endpoint","cliffs_delta_CT_minus_reference"]
    with open(a.tests_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=test_fields)
        writer.writeheader()
        writer.writerows(test_rows)

    make_box_scatter(by_pop, "FROH_GE1_5MB", "FROH (ROH >=1.5 Mb)", a.froh_png, a.froh_pdf)
    make_box_scatter(by_pop, "FROH_GE5MB", "FROH from ROH >=5 Mb", a.long_png, a.long_pdf)
    make_scatter(by_pop, a.scatter_png, a.scatter_pdf)

    with open(a.readme_out, "w", encoding="utf-8") as handle:
        handle.write(
            "ROH/FROH population distributions\n"
            "=================================\n\n"
            "Primary endpoint: FROH from all ROH >=1.5 Mb.\n"
            "Secondary endpoint: FROH contributed by ROH >=5 Mb.\n\n"
            "Five prespecified CT-vs-reference comparisons are tested separately for each endpoint.\n"
            "Tests are two-sided Mann-Whitney rank tests with tie-corrected normal approximation.\n"
            "Holm correction is applied across the five CT-vs-reference contrasts within each endpoint.\n"
            "Cliff's delta is reported as an effect size; positive values indicate larger values in CT.\n"
            "Means, medians, quartiles and individual-level plots are retained because ROH distributions\n"
            "can be strongly right-skewed and a mean alone can be driven by highly autozygous individuals.\n"
        )

if __name__ == "__main__":
    main()
