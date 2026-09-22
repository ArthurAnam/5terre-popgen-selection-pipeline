#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path


def percentile(sorted_values, p):
    if not sorted_values:
        return float("nan")
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    x = (len(sorted_values) - 1) * p
    lo = int(math.floor(x))
    hi = int(math.ceil(x))
    if lo == hi:
        return float(sorted_values[lo])
    frac = x - lo
    return sorted_values[lo] * (1.0 - frac) + sorted_values[hi] * frac


def read_next_selected(handle):
    for raw in handle:
        value = raw.strip()
        if value:
            return value
    return None


def flush_chromosome(chrom, positions, thresholds_bp, rows, global_counts):
    if chrom is None or not positions:
        return
    gaps = [b - a for a, b in zip(positions, positions[1:])]
    sorted_gaps = sorted(gaps)
    first_pos = positions[0]
    last_pos = positions[-1]
    span_bp = last_pos - first_pos + 1
    n = len(positions)
    mean_gap = (sum(gaps) / len(gaps)) if gaps else float("nan")
    median_gap = percentile(sorted_gaps, 0.5)
    p95_gap = percentile(sorted_gaps, 0.95)
    p99_gap = percentile(sorted_gaps, 0.99)
    max_gap = max(gaps) if gaps else 0

    rows.append(
        {
            "chromosome": chrom,
            "n_markers": n,
            "first_position_bp": first_pos,
            "last_position_bp": last_pos,
            "terminal_marker_span_bp": span_bp,
            "markers_per_mb_terminal_span": n / (span_bp / 1_000_000.0),
            "mean_intermarker_gap_bp": mean_gap,
            "median_intermarker_gap_bp": median_gap,
            "p95_intermarker_gap_bp": p95_gap,
            "p99_intermarker_gap_bp": p99_gap,
            "max_intermarker_gap_bp": max_gap,
            "n_intermarker_gaps": len(gaps),
        }
    )

    for gap in gaps:
        global_counts["n_gaps"] += 1
        global_counts["sum_gap_bp"] += gap
        if gap > global_counts["max_gap_bp"]:
            global_counts["max_gap_bp"] = gap
        for threshold in thresholds_bp:
            if gap > threshold:
                global_counts["threshold_counts"][threshold] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pvar", required=True)
    parser.add_argument("--snplist", required=True)
    parser.add_argument("--maf-threshold", type=float, required=True)
    parser.add_argument("--expected-variants", type=int, required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    parser.add_argument("--panel-label", default="joint_harmonized_ct_plus_1kg_eur")
    parser.add_argument("--ld-pruning-label", default="NO")
    parser.add_argument("--gap-thresholds-kb", required=True)
    parser.add_argument("--by-chromosome-out", required=True)
    parser.add_argument("--gap-thresholds-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--readme-out", required=True)
    args = parser.parse_args()

    thresholds_kb = [int(x) for x in args.gap_thresholds_kb.split(",") if x]
    thresholds_bp = [x * 1000 for x in thresholds_kb]

    rows = []
    global_counts = {
        "n_gaps": 0,
        "sum_gap_bp": 0,
        "max_gap_bp": 0,
        "threshold_counts": {x: 0 for x in thresholds_bp},
    }

    selected_count = 0
    current_chrom = None
    positions = []
    last_selected_id = None

    with open(args.snplist, "r", encoding="utf-8") as selected_handle, open(
        args.pvar, "r", encoding="utf-8"
    ) as pvar_handle:
        selected_id = read_next_selected(selected_handle)

        for raw in pvar_handle:
            if selected_id is None:
                break
            if not raw.strip() or raw.startswith("##"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if fields[0].startswith("#"):
                continue
            if len(fields) < 3:
                raise SystemExit("ERROR: malformed PVAR row")
            chrom, pos_text, variant_id = fields[0], fields[1], fields[2]

            if variant_id != selected_id:
                continue

            if variant_id == last_selected_id:
                raise SystemExit(f"ERROR: duplicate selected variant ID: {variant_id}")
            last_selected_id = variant_id

            if chrom not in {str(x) for x in range(1, 23)}:
                raise SystemExit(f"ERROR: unexpected chromosome in ROH panel: {chrom}")
            pos = int(pos_text)

            if current_chrom is None:
                current_chrom = chrom
            elif chrom != current_chrom:
                flush_chromosome(
                    current_chrom,
                    positions,
                    thresholds_bp,
                    rows,
                    global_counts,
                )
                current_chrom = chrom
                positions = []

            if positions and pos <= positions[-1]:
                raise SystemExit(
                    f"ERROR: non-increasing positions on chromosome {chrom}: "
                    f"{positions[-1]} then {pos}"
                )
            positions.append(pos)
            selected_count += 1
            selected_id = read_next_selected(selected_handle)

        flush_chromosome(
            current_chrom,
            positions,
            thresholds_bp,
            rows,
            global_counts,
        )

        if selected_id is not None:
            raise SystemExit(
                f"ERROR: selected variant {selected_id} was not found in PVAR in expected order"
            )

    if selected_count != args.expected_variants:
        raise SystemExit(
            f"ERROR: expected {args.expected_variants} MAF-filtered variants, "
            f"found {selected_count}"
        )
    if len(rows) != 22:
        raise SystemExit(f"ERROR: expected 22 autosomes, found {len(rows)}")

    outdir = Path(args.summary_out).parent
    outdir.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "chromosome",
        "n_markers",
        "first_position_bp",
        "last_position_bp",
        "terminal_marker_span_bp",
        "markers_per_mb_terminal_span",
        "mean_intermarker_gap_bp",
        "median_intermarker_gap_bp",
        "p95_intermarker_gap_bp",
        "p99_intermarker_gap_bp",
        "max_intermarker_gap_bp",
        "n_intermarker_gaps",
    ]
    with open(args.by_chromosome_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with open(args.gap_thresholds_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gap_threshold_kb", "n_gaps_exceeding", "fraction_of_all_gaps"])
        for threshold_kb, threshold_bp in zip(thresholds_kb, thresholds_bp):
            count = global_counts["threshold_counts"][threshold_bp]
            fraction = count / global_counts["n_gaps"] if global_counts["n_gaps"] else 0.0
            writer.writerow([threshold_kb, count, fraction])

    total_span_bp = sum(row["terminal_marker_span_bp"] for row in rows)
    total_markers = sum(row["n_markers"] for row in rows)
    mean_gap = (
        global_counts["sum_gap_bp"] / global_counts["n_gaps"]
        if global_counts["n_gaps"]
        else float("nan")
    )

    with open(args.summary_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["expected_joint_samples", args.expected_samples])
        writer.writerow(["maf_threshold", args.maf_threshold])
        writer.writerow(["maf_scope", args.panel_label])
        writer.writerow(["ld_pruning_before_density_audit", args.ld_pruning_label])
        writer.writerow(["retained_variants", total_markers])
        writer.writerow(["autosomes", len(rows)])
        writer.writerow(["sum_terminal_marker_spans_bp", total_span_bp])
        writer.writerow(
            ["markers_per_mb_across_terminal_marker_spans", total_markers / (total_span_bp / 1_000_000.0)]
        )
        writer.writerow(["mean_intermarker_gap_bp_genomewide", mean_gap])
        writer.writerow(["maximum_intermarker_gap_bp_genomewide", global_counts["max_gap_bp"]])
        writer.writerow(["total_intermarker_gaps", global_counts["n_gaps"]])

    with open(args.readme_out, "w", encoding="utf-8") as handle:
        handle.write(
            "ROH common-marker density audit\n"
            "===============================\n\n"
            f"Joint samples: {args.expected_samples}\n"
            f"MAF threshold: {args.maf_threshold}\n"
            f"Retained variants: {total_markers}\n"
            f"Panel label: {args.panel_label}\n"
            f"LD pruning before this audit: {args.ld_pruning_label}\n\n"
            "This step does not call ROH. It measures the physical density and gap\n"
            "structure of the common harmonized marker panel before minimum SNP\n"
            "count, minimum ROH length, density, gap, and LD-pruning parameters are\n"
            "frozen. The summed terminal-marker span is a density descriptor only;\n"
            "it is not the final denominator for FROH.\n"
        )


if __name__ == "__main__":
    main()
