#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ld-gz", required=True)
    p.add_argument("--maf-label", required=True)
    p.add_argument("--baseline-kb", type=float, required=True)
    p.add_argument("--baseline-half-width-kb", type=float, required=True)
    p.add_argument("--decay-fraction", type=float, required=True)
    p.add_argument("--bin-width-kb", type=float, required=True)
    p.add_argument("--max-distance-kb", type=float, required=True)
    p.add_argument("--stable-bins", type=int, required=True)
    p.add_argument("--bins-out", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--plot-png", required=True)
    p.add_argument("--plot-pdf", required=True)
    p.add_argument("--readme-out", required=True)
    a = p.parse_args()

    bin_width_bp = int(round(a.bin_width_kb * 1000.0))
    max_bp = int(round(a.max_distance_kb * 1000.0))
    n_bins = max_bp // bin_width_bp + 1
    counts = [0] * n_bins
    sums = [0.0] * n_bins
    sums_sq = [0.0] * n_bins

    baseline_lo = int(round((a.baseline_kb - a.baseline_half_width_kb) * 1000.0))
    baseline_hi = int(round((a.baseline_kb + a.baseline_half_width_kb) * 1000.0))
    baseline_n = 0
    baseline_sum = 0.0
    total_pairs = 0

    with gzip.open(a.ld_gz, "rt", encoding="utf-8") as handle:
        header = None
        for raw in handle:
            if not raw.strip():
                continue
            f = raw.split()
            if header is None:
                header = f
                required = {"CHR_A","BP_A","SNP_A","CHR_B","BP_B","SNP_B","R2"}
                if not required.issubset(set(header)):
                    raise SystemExit(f"ERROR: unexpected PLINK LD header: {header}")
                continue
            row = dict(zip(header, f))
            if row["CHR_A"] != row["CHR_B"] or row["SNP_A"] == row["SNP_B"]:
                continue
            dist = abs(int(row["BP_B"]) - int(row["BP_A"]))
            if dist <= 0 or dist > max_bp:
                continue
            r2 = float(row["R2"])
            if not math.isfinite(r2):
                continue
            idx = dist // bin_width_bp
            if idx >= n_bins:
                continue
            counts[idx] += 1
            sums[idx] += r2
            sums_sq[idx] += r2 * r2
            total_pairs += 1
            if baseline_lo <= dist < baseline_hi:
                baseline_n += 1
                baseline_sum += r2

    if baseline_n == 0:
        raise SystemExit("ERROR: no LD pairs in the baseline distance interval")
    baseline_mean = baseline_sum / baseline_n
    threshold = baseline_mean * a.decay_fraction

    rows = []
    means = []
    for i in range(n_bins):
        if counts[i]:
            m = sums[i] / counts[i]
            var = max(0.0, sums_sq[i] / counts[i] - m * m)
            sd = math.sqrt(var)
            se = sd / math.sqrt(counts[i])
        else:
            m = float("nan")
            sd = float("nan")
            se = float("nan")
        lower = i * bin_width_bp
        upper = min(max_bp, (i + 1) * bin_width_bp)
        center = (lower + upper) / 2.0
        means.append(m)
        rows.append({
            "bin_lower_kb": lower / 1000.0,
            "bin_upper_kb": upper / 1000.0,
            "bin_center_kb": center / 1000.0,
            "n_pairs": counts[i],
            "mean_r2": m,
            "sd_r2": sd,
            "se_r2": se,
        })

    search_start_bp = baseline_hi
    first_crossing = None
    stable_crossing = None
    for i, row in enumerate(rows):
        if row["bin_lower_kb"] * 1000.0 < search_start_bp:
            continue
        if row["n_pairs"] > 0 and math.isfinite(row["mean_r2"]) and row["mean_r2"] < threshold:
            if first_crossing is None:
                first_crossing = row["bin_center_kb"]
            ok = True
            for j in range(i, min(i + a.stable_bins, len(rows))):
                if rows[j]["n_pairs"] == 0 or not math.isfinite(rows[j]["mean_r2"]) or rows[j]["mean_r2"] >= threshold:
                    ok = False
                    break
            if ok and i + a.stable_bins <= len(rows):
                stable_crossing = row["bin_center_kb"]
                break

    with open(a.bins_out, "w", encoding="utf-8", newline="") as out:
        fields = ["bin_lower_kb","bin_upper_kb","bin_center_kb","n_pairs","mean_r2","sd_r2","se_r2"]
        w = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    with open(a.summary_out, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["metric","value"])
        w.writerow(["maf_label", a.maf_label])
        w.writerow(["total_valid_ld_pairs", total_pairs])
        w.writerow(["baseline_distance_kb", a.baseline_kb])
        w.writerow(["baseline_interval_kb", f"{a.baseline_kb-a.baseline_half_width_kb:g}-{a.baseline_kb+a.baseline_half_width_kb:g}"])
        w.writerow(["baseline_pair_count", baseline_n])
        w.writerow(["baseline_mean_r2", baseline_mean])
        w.writerow(["decay_fraction", a.decay_fraction])
        w.writerow(["decay_threshold_r2", threshold])
        w.writerow(["first_crossing_bin_center_kb", "NA" if first_crossing is None else first_crossing])
        w.writerow(["stable_crossing_bin_center_kb", "NA" if stable_crossing is None else stable_crossing])
        w.writerow(["stable_crossing_required_consecutive_bins", a.stable_bins])
        w.writerow(["max_distance_kb", a.max_distance_kb])

    xs = [r["bin_center_kb"] for r in rows if r["n_pairs"] > 0]
    ys = [r["mean_r2"] for r in rows if r["n_pairs"] > 0]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xs, ys, linewidth=1.2)
    ax.axhline(threshold, linestyle="--", linewidth=1.0)
    if first_crossing is not None:
        ax.axvline(first_crossing, linestyle=":", linewidth=1.0)
    ax.set_xlabel("Physical distance between SNPs (kb)")
    ax.set_ylabel("Mean pairwise r²")
    ax.set_xlim(0, a.max_distance_kb)
    ax.set_ylim(bottom=0)
    ax.set_title(f"Cinque Terre LD decay ({a.maf_label})")
    fig.tight_layout()
    fig.savefig(a.plot_png, dpi=300)
    fig.savefig(a.plot_pdf)
    plt.close(fig)

    with open(a.readme_out, "w", encoding="utf-8") as out:
        out.write(
            "Selection LD-decay audit\n"
            "========================\n\n"
            f"Panel: final QCed Cinque Terre, {a.maf_label}.\n"
            "Statistic: PLINK 1.9 unphased hard-call r^2.\n"
            f"Pairs are evaluated up to {a.max_distance_kb:g} kb from deterministic genome-wide anchors.\n"
            f"Baseline r^2: mean among pairs separated by {a.baseline_kb-a.baseline_half_width_kb:g}-"
            f"{a.baseline_kb+a.baseline_half_width_kb:g} kb.\n"
            f"Decay threshold: {a.decay_fraction:g} x baseline mean r^2.\n"
            "The first crossing is the direct literature-matched physical-decay estimate.\n"
            f"A {a.stable_bins}-consecutive-bin crossing is also reported only as a stability diagnostic.\n"
            "No LASSI SNP-window size is frozen at this stage.\n"
        )

if __name__ == "__main__":
    main()
