#!/usr/bin/env python3
import argparse
import bisect
import csv
from collections import defaultdict
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bim", required=True)
    p.add_argument("--anchors-per-chromosome", type=int, required=True)
    p.add_argument("--edge-margin-kb", type=float, required=True)
    p.add_argument("--anchors-out", required=True)
    p.add_argument("--summary-out", required=True)
    args = p.parse_args()

    by_chr = defaultdict(list)
    with open(args.bim, "r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            f = raw.split()
            if len(f) < 4:
                raise SystemExit("ERROR: malformed BIM row")
            chrom, vid, pos = f[0], f[1], int(f[3])
            if chrom not in {str(i) for i in range(1, 23)}:
                continue
            by_chr[chrom].append((pos, vid))

    margin = int(round(args.edge_margin_kb * 1000.0))
    selected_all = []
    summary = []

    for chrom in [str(i) for i in range(1, 23)]:
        rows = sorted(by_chr.get(chrom, []))
        if not rows:
            raise SystemExit(f"ERROR: chromosome {chrom} has no variants")
        min_pos, max_pos = rows[0][0], rows[-1][0]
        lo, hi = min_pos + margin, max_pos - margin
        eligible = [(pos, vid) for pos, vid in rows if lo <= pos <= hi]
        if not eligible:
            raise SystemExit(f"ERROR: chromosome {chrom} has no variants after edge margin")

        positions = [x[0] for x in eligible]
        n_target = min(args.anchors_per_chromosome, len(eligible))
        if n_target == 1:
            targets = [(lo + hi) // 2]
        else:
            targets = [round(lo + i * (hi - lo) / (n_target - 1)) for i in range(n_target)]

        chosen = []
        used = set()
        for target in targets:
            idx = bisect.bisect_left(positions, target)
            candidates = []
            if idx < len(positions):
                candidates.append(idx)
            if idx > 0:
                candidates.append(idx - 1)
            candidates.sort(key=lambda j: (abs(positions[j] - target), positions[j]))
            pick = next((j for j in candidates if eligible[j][1] not in used), None)
            if pick is None:
                # Dense WGS data make this unlikely; search locally if necessary.
                radius = 1
                while pick is None and (idx - radius >= 0 or idx + radius < len(eligible)):
                    for j in (idx - radius, idx + radius):
                        if 0 <= j < len(eligible) and eligible[j][1] not in used:
                            pick = j
                            break
                    radius += 1
            if pick is None:
                continue
            used.add(eligible[pick][1])
            chosen.append(eligible[pick])

        chosen.sort()
        selected_all.extend((chrom, pos, vid) for pos, vid in chosen)
        spacings = [chosen[i][0] - chosen[i-1][0] for i in range(1, len(chosen))]
        summary.append({
            "chromosome": chrom,
            "total_variants": len(rows),
            "eligible_variants": len(eligible),
            "selected_anchors": len(chosen),
            "first_anchor_bp": chosen[0][0],
            "last_anchor_bp": chosen[-1][0],
            "mean_anchor_spacing_bp": (sum(spacings) / len(spacings)) if spacings else 0.0,
        })

    Path(args.anchors_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.anchors_out, "w", encoding="utf-8") as out:
        for chrom, pos, vid in selected_all:
            out.write(vid + "\n")

    with open(args.summary_out, "w", encoding="utf-8", newline="") as out:
        fields = ["chromosome","total_variants","eligible_variants","selected_anchors",
                  "first_anchor_bp","last_anchor_bp","mean_anchor_spacing_bp"]
        w = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(summary)

if __name__ == "__main__":
    main()
