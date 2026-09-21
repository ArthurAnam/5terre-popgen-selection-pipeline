#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

DUPLICATE_MIN = 0.354
FIRST_MIN = 0.177
SECOND_MIN = 0.0884
THIRD_MIN = 0.0442


def parse_args():
    p = argparse.ArgumentParser(description="Summarize and plot KING pairwise kinship results.")
    p.add_argument("--kin0", required=True)
    p.add_argument("--fam", required=True)
    p.add_argument("--pairs-out", required=True)
    p.add_argument("--candidates-out", required=True)
    p.add_argument("--counts-out", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--individual-summary-out", required=True)
    p.add_argument("--hist-out", required=True)
    p.add_argument("--scatter-out", required=True)
    p.add_argument("--scatter-pdf-out", required=True)
    p.add_argument("--heatmap-out", required=True)
    return p.parse_args()


def classify(k):
    if k > DUPLICATE_MIN:
        return "duplicate_or_MZ"
    if k >= FIRST_MIN:
        return "first_degree"
    if k >= SECOND_MIN:
        return "second_degree"
    if k >= THIRD_MIN:
        return "possible_third_degree_candidate"
    return "unrelated_or_more_distant"


def read_samples(path):
    samples = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError("Malformed PLINK .fam line")
            samples.append(fields[1])
    if len(samples) != len(set(samples)):
        raise ValueError("Duplicate sample IDs in PLINK .fam")
    return samples


def read_kin0(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        header = None
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if header is None:
                header = fields
                required = {"ID1", "ID2", "N_SNP", "HetHet", "IBS0", "Kinship"}
                if not required.issubset(set(header)):
                    missing = sorted(required - set(header))
                    raise ValueError(f"KING .kin0 missing required columns: {missing}")
                idx = {name: header.index(name) for name in required}
                continue
            rows.append({
                "sample1": fields[idx["ID1"]],
                "sample2": fields[idx["ID2"]],
                "n_snp": int(fields[idx["N_SNP"]]),
                "hethet": float(fields[idx["HetHet"]]),
                "ibs0": float(fields[idx["IBS0"]]),
                "kinship": float(fields[idx["Kinship"]]),
            })
    if header is None:
        raise ValueError("Empty KING .kin0 file")
    return rows


def write_tables(rows, pairs_out, candidates_out, counts_out, summary_out, n_samples):
    for row in rows:
        row["relationship"] = classify(row["kinship"])

    rows_sorted = sorted(rows, key=lambda r: r["kinship"], reverse=True)
    fields = ["sample1", "sample2", "n_snp", "hethet", "ibs0", "kinship", "relationship"]

    with open(pairs_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_sorted)

    candidates = [r for r in rows_sorted if r["kinship"] >= THIRD_MIN]
    robust_close = [r for r in rows_sorted if r["kinship"] >= SECOND_MIN]
    with open(candidates_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(candidates)

    order = [
        "duplicate_or_MZ",
        "first_degree",
        "second_degree",
        "possible_third_degree_candidate",
        "unrelated_or_more_distant",
    ]
    counts = Counter(r["relationship"] for r in rows)
    with open(counts_out, "w", encoding="utf-8") as handle:
        handle.write("relationship_class\tpair_count\n")
        for key in order:
            handle.write(f"{key}\t{counts[key]}\n")

    expected_pairs = n_samples * (n_samples - 1) // 2
    if len(rows) != expected_pairs:
        raise ValueError(
            f"Expected {expected_pairs} pairwise comparisons for {n_samples} samples, "
            f"observed {len(rows)}"
        )

    kinships = np.array([r["kinship"] for r in rows], dtype=float)
    with open(summary_out, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"n_samples\t{n_samples}\n")
        handle.write(f"n_pairwise_comparisons\t{len(rows)}\n")
        handle.write(f"robust_candidate_pairs_second_degree_or_closer\t{len(robust_close)}\n")
        handle.write(f"exploratory_pairs_at_third_degree_threshold_or_closer\t{len(candidates)}\n")
        handle.write(f"max_kinship\t{np.max(kinships):.8f}\n")
        handle.write(f"min_kinship\t{np.min(kinships):.8f}\n")
        handle.write(f"mean_kinship\t{np.mean(kinships):.8f}\n")
        handle.write(f"median_kinship\t{np.median(kinships):.8f}\n")
        handle.write(f"second_degree_threshold\t{SECOND_MIN}\n")
        handle.write(f"exploratory_third_degree_threshold\t{THIRD_MIN}\n")
        handle.write("automatic_sample_exclusion\tNO\n")


def write_individual_summary(rows, samples, out_path):
    values = defaultdict(list)
    for row in rows:
        values[row["sample1"]].append(row["kinship"])
        values[row["sample2"]].append(row["kinship"])

    expected_n = len(samples) - 1
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("sample\tn_pairs\tmean_kinship\tmedian_kinship\tmin_kinship\tmax_kinship\n")
        for sample in samples:
            x = np.array(values[sample], dtype=float)
            if len(x) != expected_n:
                raise ValueError(
                    f"Sample {sample}: expected {expected_n} pairwise values, observed {len(x)}"
                )
            handle.write(
                f"{sample}\t{len(x)}\t{np.mean(x):.8f}\t{np.median(x):.8f}\t"
                f"{np.min(x):.8f}\t{np.max(x):.8f}\n"
            )


def make_plots(rows, samples, hist_out, scatter_out, scatter_pdf_out, heatmap_out):
    kinship = np.array([r["kinship"] for r in rows], dtype=float)
    ibs0 = np.array([r["ibs0"] for r in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(kinship, bins=50)
    ax.axvline(0, linewidth=1)
    for threshold in [THIRD_MIN, SECOND_MIN, FIRST_MIN, DUPLICATE_MIN]:
        ax.axvline(threshold, linestyle="--", linewidth=1)
    ax.set_xlabel("KING-Robust kinship estimate")
    ax.set_ylabel("Pair count")
    ax.set_title("Cinque Terre pairwise KING-Robust kinship")
    fig.tight_layout()
    fig.savefig(hist_out, dpi=180)
    plt.close(fig)

    # Publication-oriented KING diagnostic.  Only the third- and second-degree
    # boundaries are drawn because all observed pairs are far below the
    # first-degree and duplicate/MZ ranges.  Showing the much higher boundaries
    # would compress the observed data into a narrow strip at the bottom.
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    ax.scatter(ibs0, kinship, s=20, alpha=0.55, linewidths=0)

    x_span = float(np.max(ibs0) - np.min(ibs0))
    x_pad = 0.04 * x_span if x_span > 0 else 0.001
    x_min = float(np.min(ibs0) - x_pad)
    x_max = float(np.max(ibs0) + x_pad)
    ax.set_xlim(x_min, x_max)

    y_min = min(-0.06, float(np.min(kinship) - 0.005))
    y_max = max(0.105, SECOND_MIN + 0.012)
    ax.set_ylim(y_min, y_max)

    ax.axhline(0, linewidth=0.9, color="0.35")
    for threshold in [THIRD_MIN, SECOND_MIN]:
        ax.axhline(threshold, linestyle="--", linewidth=1.0, color="0.35")

    label_x = x_max - 0.015 * (x_max - x_min)
    ax.text(
        label_x, THIRD_MIN + 0.002,
        "third-degree screening threshold (0.0442)",
        ha="right", va="bottom", fontsize=9,
    )
    ax.text(
        label_x, SECOND_MIN + 0.002,
        "second-degree lower bound (0.0884)",
        ha="right", va="bottom", fontsize=9,
    )

    max_idx = int(np.argmax(kinship))
    ax.scatter(
        [ibs0[max_idx]], [kinship[max_idx]],
        s=58, facecolors="none", edgecolors="black", linewidths=1.0, zorder=3,
    )
    ax.annotate(
        f"max observed = {kinship[max_idx]:.4f}",
        xy=(ibs0[max_idx], kinship[max_idx]),
        xytext=(12, 12),
        textcoords="offset points",
        fontsize=9,
        arrowprops={"arrowstyle": "-", "linewidth": 0.8},
    )

    ax.set_xlabel("IBS0 proportion")
    ax.set_ylabel("KING-Robust kinship estimate")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(scatter_out, dpi=300, bbox_inches="tight")
    fig.savefig(scatter_pdf_out, bbox_inches="tight")
    plt.close(fig)

    n = len(samples)
    sample_index = {sample: i for i, sample in enumerate(samples)}
    matrix = np.full((n, n), np.nan, dtype=float)

    for row in rows:
        i = sample_index[row["sample1"]]
        j = sample_index[row["sample2"]]
        matrix[i, j] = row["kinship"]
        matrix[j, i] = row["kinship"]

    pair_min = float(np.nanmin(matrix))
    pair_max = float(np.nanmax(matrix))
    norm = TwoSlopeNorm(vmin=pair_min, vcenter=0, vmax=pair_max) if pair_min < 0 < pair_max else None
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("white")

    fig, ax = plt.subplots(figsize=(13, 11))
    image = ax.imshow(np.ma.masked_invalid(matrix), aspect="auto", cmap=cmap, norm=norm)
    ax.set_title("Pairwise KING-Robust kinship matrix (diagonal omitted)")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(samples, rotation=90, fontsize=6)
    ax.set_yticklabels(samples, fontsize=6)
    fig.colorbar(image, ax=ax, label="KING-Robust kinship estimate")
    fig.tight_layout()
    fig.savefig(heatmap_out, dpi=180)
    plt.close(fig)


def main():
    args = parse_args()
    samples = read_samples(args.fam)
    rows = read_kin0(args.kin0)
    write_tables(
        rows,
        args.pairs_out,
        args.candidates_out,
        args.counts_out,
        args.summary_out,
        len(samples),
    )
    write_individual_summary(rows, samples, args.individual_summary_out)
    make_plots(
        rows,
        samples,
        args.hist_out,
        args.scatter_out,
        args.scatter_pdf_out,
        args.heatmap_out,
    )


if __name__ == "__main__":
    main()
