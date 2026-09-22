#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt


POP_ORDER = ["CEU", "FIN", "GBR", "IBS", "TSI", "CT"]
MARKERS = {
    "CEU": "o",
    "FIN": "o",
    "GBR": "o",
    "IBS": "o",
    "TSI": "o",
    "CT": "*",
}


def read_annotations(path: Path):
    mapping = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            mapping[row["sample_id"]] = row["population"]
    return mapping


def read_eval(path: Path):
    values = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line:
                values.append(float(line))
    if not values:
        raise SystemExit("ERROR: empty smartpca eigenvalue file")
    positive_sum = sum(v for v in values if v > 0)
    if positive_sum <= 0:
        raise SystemExit("ERROR: non-positive sum of smartpca eigenvalues")
    pve = [100.0 * v / positive_sum for v in values]
    return values, pve


def read_evec(path: Path):
    rows = []
    n_pcs = None
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise SystemExit(f"ERROR: malformed smartpca evec row: {line}")
            current_n = len(fields) - 2
            if n_pcs is None:
                n_pcs = current_n
            elif current_n != n_pcs:
                raise SystemExit("ERROR: inconsistent PC count in evec")
            sample = fields[0]
            pcs = [float(x) for x in fields[1:-1]]
            pop = fields[-1]
            rows.append((sample, pcs, pop))
    if not rows:
        raise SystemExit("ERROR: empty smartpca evec file")
    return rows, n_pcs


def make_plot(rows, pve, x_index, y_index, png_path):
    fig, ax = plt.subplots(figsize=(7.0, 5.8))

    for pop in POP_ORDER:
        subset = [(pcs[x_index], pcs[y_index]) for _, pcs, p in rows if p == pop]
        if not subset:
            continue
        xs = [x for x, _ in subset]
        ys = [y for _, y in subset]
        if pop == "CT":
            ax.scatter(
                xs, ys, label="Cinque Terre", marker=MARKERS[pop],
                s=70, linewidths=0.7, edgecolors="black", zorder=5
            )
        else:
            ax.scatter(
                xs, ys, label=pop, marker=MARKERS[pop],
                s=24, alpha=0.72, linewidths=0
            )

    ax.axhline(0, linewidth=0.5, alpha=0.35)
    ax.axvline(0, linewidth=0.5, alpha=0.35)
    ax.set_xlabel(f"PC{x_index + 1} ({pve[x_index]:.2f}%)")
    ax.set_ylabel(f"PC{y_index + 1} ({pve[y_index]:.2f}%)")
    ax.legend(frameon=False, fontsize=9, ncol=2)
    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evec", required=True)
    parser.add_argument("--eval", required=True)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    parser.add_argument("--coordinates-out", required=True)
    parser.add_argument("--eigenvalues-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--pc12-png", required=True)
    parser.add_argument("--pc23-png", required=True)
    args = parser.parse_args()

    annotations = read_annotations(Path(args.annotations))
    evals, pve = read_eval(Path(args.eval))
    rows, n_pcs = read_evec(Path(args.evec))

    if len(rows) != args.expected_samples:
        raise SystemExit(
            f"ERROR: expected {args.expected_samples} smartpca samples, found {len(rows)}"
        )
    if n_pcs < 3:
        raise SystemExit("ERROR: at least three PCs are required for the requested plots")

    seen = set()
    for sample, _, pop in rows:
        if sample in seen:
            raise SystemExit(f"ERROR: duplicate smartpca sample {sample}")
        seen.add(sample)
        expected_pop = annotations.get(sample)
        if expected_pop is None:
            raise SystemExit(f"ERROR: smartpca sample {sample} missing from annotations")
        if pop != expected_pop:
            raise SystemExit(
                f"ERROR: population label mismatch for {sample}: {pop} != {expected_pop}"
            )

    if set(annotations) != seen:
        raise SystemExit("ERROR: smartpca sample set differs from annotation sample set")

    coordinates_out = Path(args.coordinates_out)
    eigenvalues_out = Path(args.eigenvalues_out)
    summary_out = Path(args.summary_out)
    coordinates_out.parent.mkdir(parents=True, exist_ok=True)

    with coordinates_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "population"] + [f"PC{i}" for i in range(1, n_pcs + 1)])
        for sample, pcs, pop in rows:
            writer.writerow([sample, pop] + pcs)

    with eigenvalues_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["PC", "eigenvalue", "variance_percent"])
        for i, (eigenvalue, variance) in enumerate(zip(evals, pve), start=1):
            writer.writerow([i, eigenvalue, variance])

    counts = Counter(pop for _, _, pop in rows)
    with summary_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["panel", args.panel])
        writer.writerow(["samples", len(rows)])
        writer.writerow(["pcs_output", n_pcs])
        for pop in POP_ORDER:
            writer.writerow([f"samples_{pop}", counts[pop]])
        for i in range(min(n_pcs, 10)):
            writer.writerow([f"PC{i+1}_eigenvalue", evals[i]])
            writer.writerow([f"PC{i+1}_variance_percent", pve[i]])

    make_plot(rows, pve, 0, 1, args.pc12_png)
    make_plot(rows, pve, 1, 2, args.pc23_png)


if __name__ == "__main__":
    main()
