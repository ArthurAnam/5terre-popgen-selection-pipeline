#!/usr/bin/env python3
import argparse
import bisect
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

def read_bim(path):
    by_chr = defaultdict(list)
    id_to_chr_pos = {}
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            f = raw.split()
            if len(f) < 4:
                raise SystemExit("ERROR: malformed BIM row")
            chrom, vid, pos = f[0], f[1], int(f[3])
            if chrom not in {str(i) for i in range(1, 23)}:
                continue
            by_chr[chrom].append(pos)
            id_to_chr_pos[vid] = (chrom, pos)
    for chrom in by_chr:
        by_chr[chrom].sort()
    return by_chr, id_to_chr_pos

def read_anchors(path, id_to_chr_pos):
    anchors = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            vid = raw.strip()
            if not vid:
                continue
            if vid not in id_to_chr_pos:
                raise SystemExit(f"ERROR: anchor {vid} not found in BIM")
            chrom, pos = id_to_chr_pos[vid]
            anchors.append((vid, chrom, pos))
    if not anchors:
        raise SystemExit("ERROR: no anchors")
    return anchors

def count_in_window(positions, center, width_bp):
    half = width_bp / 2.0
    lo = center - half
    hi = center + half
    left = bisect.bisect_left(positions, lo)
    right = bisect.bisect_right(positions, hi)
    return right - left

def q(values, frac):
    vals = sorted(values)
    if len(vals) == 1:
        return float(vals[0])
    pos = (len(vals) - 1) * frac
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(vals[lo])
    w = pos - lo
    return vals[lo] * (1.0 - w) + vals[hi] * w

def round_half_up(x):
    return int(math.floor(x + 0.5))

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bim", required=True)
    p.add_argument("--anchors", required=True)
    p.add_argument("--anchor-bim", default=None, help="Optional BIM used only to resolve anchor IDs/coordinates; defaults to --bim.")
    p.add_argument("--primary-width-kb", type=float, required=True)
    p.add_argument("--stable-width-kb", type=float, required=True)
    p.add_argument("--step-fraction", type=float, required=True)
    p.add_argument("--anchor-counts-out", required=True)
    p.add_argument("--chrom-summary-out", required=True)
    p.add_argument("--summary-out", required=True)
    args = p.parse_args()

    by_chr, id_map = read_bim(args.bim)
    if args.anchor_bim:
        _, anchor_id_map = read_bim(args.anchor_bim)
    else:
        anchor_id_map = id_map
    anchors = read_anchors(args.anchors, anchor_id_map)

    primary_bp = args.primary_width_kb * 1000.0
    stable_bp = args.stable_width_kb * 1000.0

    rows = []
    for vid, chrom, pos in anchors:
        primary_n = count_in_window(by_chr[chrom], pos, primary_bp)
        stable_n = count_in_window(by_chr[chrom], pos, stable_bp)
        rows.append({
            "anchor_id": vid,
            "chromosome": chrom,
            "position_bp": pos,
            "primary_width_kb": args.primary_width_kb,
            "primary_snp_count": primary_n,
            "stable_width_kb": args.stable_width_kb,
            "stable_snp_count": stable_n,
        })

    Path(args.anchor_counts_out).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "anchor_id","chromosome","position_bp",
        "primary_width_kb","primary_snp_count",
        "stable_width_kb","stable_snp_count"
    ]
    with open(args.anchor_counts_out, "w", encoding="utf-8", newline="") as out:
        w = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    with open(args.chrom_summary_out, "w", encoding="utf-8", newline="") as out:
        fields2 = [
            "chromosome","n_anchors",
            "primary_mean_snps","primary_median_snps","primary_q1_snps","primary_q3_snps",
            "stable_mean_snps","stable_median_snps","stable_q1_snps","stable_q3_snps"
        ]
        w = csv.DictWriter(out, delimiter="\t", fieldnames=fields2)
        w.writeheader()
        for chrom in [str(i) for i in range(1, 23)]:
            sub = [r for r in rows if r["chromosome"] == chrom]
            pvals = [r["primary_snp_count"] for r in sub]
            svals = [r["stable_snp_count"] for r in sub]
            w.writerow({
                "chromosome": chrom,
                "n_anchors": len(sub),
                "primary_mean_snps": statistics.mean(pvals),
                "primary_median_snps": statistics.median(pvals),
                "primary_q1_snps": q(pvals, 0.25),
                "primary_q3_snps": q(pvals, 0.75),
                "stable_mean_snps": statistics.mean(svals),
                "stable_median_snps": statistics.median(svals),
                "stable_q1_snps": q(svals, 0.25),
                "stable_q3_snps": q(svals, 0.75),
            })

    pvals = [r["primary_snp_count"] for r in rows]
    svals = [r["stable_snp_count"] for r in rows]
    pmed = statistics.median(pvals)
    smed = statistics.median(svals)
    winsize = round_half_up(pmed)
    stable_winsize = round_half_up(smed)
    winstep = max(1, round_half_up(winsize * args.step_fraction))
    stable_winstep = max(1, round_half_up(stable_winsize * args.step_fraction))

    with open(args.summary_out, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["metric","value"])
        w.writerow(["n_anchors", len(rows)])
        w.writerow(["primary_physical_width_kb", args.primary_width_kb])
        w.writerow(["primary_mean_snp_count", statistics.mean(pvals)])
        w.writerow(["primary_median_snp_count", pmed])
        w.writerow(["primary_q1_snp_count", q(pvals, 0.25)])
        w.writerow(["primary_q3_snp_count", q(pvals, 0.75)])
        w.writerow(["primary_min_snp_count", min(pvals)])
        w.writerow(["primary_max_snp_count", max(pvals)])
        w.writerow(["lassi_winsize_snps_candidate", winsize])
        w.writerow(["lassi_winstep_snps_candidate", winstep])
        w.writerow(["stable_physical_width_kb", args.stable_width_kb])
        w.writerow(["stable_mean_snp_count", statistics.mean(svals)])
        w.writerow(["stable_median_snp_count", smed])
        w.writerow(["stable_q1_snp_count", q(svals, 0.25)])
        w.writerow(["stable_q3_snp_count", q(svals, 0.75)])
        w.writerow(["stable_min_snp_count", min(svals)])
        w.writerow(["stable_max_snp_count", max(svals)])
        w.writerow(["stable_winsize_snps_candidate", stable_winsize])
        w.writerow(["stable_winstep_snps_candidate", stable_winstep])
        w.writerow(["step_fraction", args.step_fraction])
        w.writerow(["decision_rule", "use rounded genome-wide median SNP count in the population-specific primary-width MAF>=0.05 windows; stable-width result is the stability diagnostic"])

if __name__ == "__main__":
    main()
