#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


def read_fst(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("ERROR: empty FST summary")
        normalized = [x.lstrip("#") for x in reader.fieldnames]
        rename = {old: new for old, new in zip(reader.fieldnames, normalized)}
        if "POP1" not in normalized or "POP2" not in normalized or "HUDSON_FST" not in normalized:
            raise SystemExit(f"ERROR: unexpected FST columns: {reader.fieldnames}")
        for raw in reader:
            row = {rename[k]: v for k, v in raw.items()}
            fst = float(row["HUDSON_FST"])
            if not math.isfinite(fst):
                raise SystemExit(f"ERROR: non-finite FST for {row['POP1']} vs {row['POP2']}")
            rows.append((row["POP1"], row["POP2"], fst, row.get("OBS_CT", "")))
    return rows


def average_linkage(populations, dist):
    n = len(populations)
    clusters = {i: [p] for i, p in enumerate(populations)}
    heights = {i: 0.0 for i in range(n)}
    x_positions = {i: float(i) for i in range(n)}
    leaves = {i: [i] for i in range(n)}
    merges = []
    next_id = n

    while len(clusters) > 1:
        ids = sorted(clusters)
        best = None
        for ai, a in enumerate(ids):
            for b in ids[ai + 1:]:
                vals = [dist[x][y] for x in clusters[a] for y in clusters[b]]
                d = sum(vals) / len(vals)
                key = (d, tuple(sorted(clusters[a])), tuple(sorted(clusters[b])), a, b)
                if best is None or key < best:
                    best = key
        d, _, _, a, b = best
        new_members = sorted(clusters[a] + clusters[b])
        merges.append({
            "merge_step": len(merges) + 1,
            "left_members": ",".join(sorted(clusters[a])),
            "right_members": ",".join(sorted(clusters[b])),
            "height": d,
            "new_members": ",".join(new_members),
            "left_id": a,
            "right_id": b,
            "new_id": next_id,
        })
        clusters[next_id] = new_members
        leaves[next_id] = leaves[a] + leaves[b]
        x_positions[next_id] = sum(x_positions[i] for i in leaves[next_id]) / len(leaves[next_id])
        heights[next_id] = d
        del clusters[a]
        del clusters[b]
        next_id += 1
    return merges, x_positions, heights, leaves


def draw_dendrogram(populations, merges, out):
    n = len(populations)
    x = {i: float(i) for i in range(n)}
    h = {i: 0.0 for i in range(n)}
    leaves = {i: [i] for i in range(n)}
    next_id = n

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for m in merges:
        a = m["left_id"]
        b = m["right_id"]
        height = m["height"]
        xa = sum(x[i] for i in leaves[a]) / len(leaves[a])
        xb = sum(x[i] for i in leaves[b]) / len(leaves[b])
        ax.plot([xa, xa], [h[a], height], linewidth=1.1)
        ax.plot([xb, xb], [h[b], height], linewidth=1.1)
        ax.plot([xa, xb], [height, height], linewidth=1.1)
        leaves[next_id] = leaves[a] + leaves[b]
        x[next_id] = sum(x[i] for i in leaves[next_id]) / len(leaves[next_id])
        h[next_id] = height
        next_id += 1

    ax.set_xticks(range(n))
    ax.set_xticklabels(populations)
    ax.set_ylabel("Hudson FST distance")
    ax.set_xlabel("Population")
    ax.set_title("Exploratory average-linkage clustering")
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fst-summary", required=True)
    parser.add_argument("--populations", required=True)
    parser.add_argument("--matrix-out", required=True)
    parser.add_argument("--ranked-out", required=True)
    parser.add_argument("--merges-out", required=True)
    parser.add_argument("--dendrogram-out", required=True)
    parser.add_argument("--readme-out", required=True)
    args = parser.parse_args()

    populations = [x for x in args.populations.split(",") if x]
    rows = read_fst(args.fst_summary)

    raw = {p: {q: 0.0 for q in populations} for p in populations}
    seen = set()
    for a, b, fst, _ in rows:
        if a not in raw or b not in raw:
            raise SystemExit(f"ERROR: unexpected population pair {a}, {b}")
        key = tuple(sorted((a, b)))
        if key in seen:
            raise SystemExit(f"ERROR: duplicate pair {key}")
        seen.add(key)
        raw[a][b] = fst
        raw[b][a] = fst

    expected = len(populations) * (len(populations) - 1) // 2
    if len(seen) != expected:
        raise SystemExit(f"ERROR: expected {expected} pairs, found {len(seen)}")

    # Clustering distance cannot be negative; preserve raw FST separately.
    dist = {
        p: {q: (0.0 if p == q else max(0.0, raw[p][q])) for q in populations}
        for p in populations
    }

    outdir = Path(args.matrix_out).parent
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.matrix_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["population"] + populations)
        for p in populations:
            writer.writerow([p] + [raw[p][q] for q in populations])

    ranked = sorted(rows, key=lambda x: (x[2], x[0], x[1]))
    with open(args.ranked_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "POP1", "POP2", "HUDSON_FST", "OBS_CT"])
        for i, (a, b, fst, obs) in enumerate(ranked, start=1):
            writer.writerow([i, a, b, fst, obs])

    merges, _, _, _ = average_linkage(populations, dist)
    with open(args.merges_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "merge_step",
                "left_members",
                "right_members",
                "height",
                "new_members",
                "left_id",
                "right_id",
                "new_id",
            ],
        )
        writer.writeheader()
        writer.writerows(merges)

    draw_dendrogram(populations, merges, args.dendrogram_out)

    with open(args.readme_out, "w", encoding="utf-8") as handle:
        handle.write(
            "Exploratory population grouping diagnostic\n"
            "=========================================\n\n"
            "Input: high-LD-masked, MAF-filtered, LD-pruned joint PCA panel.\n"
            "Differentiation: pairwise Hudson FST from PLINK2.\n"
            "Clustering: average linkage on pairwise FST distances.\n"
            "Negative FST estimates, if present, are retained in the matrix and\n"
            "ranked table but truncated to zero only when used as a clustering\n"
            "distance.\n\n"
            "This branch is exploratory. It does not merge population labels or\n"
            "define biological populations. If a proposed pooling becomes\n"
            "important for downstream inference, add uncertainty/stability checks\n"
            "before using it as a formal grouping.\n"
        )


if __name__ == "__main__":
    main()
