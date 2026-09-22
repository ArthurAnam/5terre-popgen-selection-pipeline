#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


POP_ORDER = ["CT", "CEU", "FIN", "GBR", "IBS", "TSI"]
REF_POPS = ["CEU", "FIN", "GBR", "IBS", "TSI"]


def read_coordinates(path: Path, n_components: int):
    data = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "population"} | {
            f"PC{i}" for i in range(1, n_components + 1)
        }
        if not required.issubset(reader.fieldnames or []):
            missing = sorted(required - set(reader.fieldnames or []))
            raise SystemExit(f"ERROR: missing coordinate columns in {path}: {missing}")
        for row in reader:
            sid = row["sample_id"]
            if sid in data:
                raise SystemExit(f"ERROR: duplicate sample ID in {path}: {sid}")
            data[sid] = {
                "population": row["population"],
                "pcs": [float(row[f"PC{i}"]) for i in range(1, n_components + 1)],
            }
    return data


def read_eigenvalues(path: Path, n_components: int):
    values = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pc = int(row["PC"])
            if pc <= n_components:
                values[pc] = {
                    "eigenvalue": float(row["eigenvalue"]),
                    "variance_percent": float(row["variance_percent"]),
                }
    if len(values) < n_components:
        raise SystemExit(
            f"ERROR: expected at least {n_components} eigenvalues in {path}, found {len(values)}"
        )
    return values


def mean(xs):
    return sum(xs) / len(xs)


def pearson(x, y):
    if len(x) != len(y) or not x:
        raise ValueError("Pearson correlation requires equal non-empty vectors")
    mx = mean(x)
    my = mean(y)
    sx = sum((v - mx) ** 2 for v in x)
    sy = sum((v - my) ** 2 for v in y)
    if sx <= 0 or sy <= 0:
        raise ValueError("Cannot correlate a constant vector")
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return cov / math.sqrt(sx * sy)


def rmse(x, y):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(x, y)) / len(x))


def centroid(data, samples, pc_indices):
    return [
        mean([data[sid]["pcs"][idx] for sid in samples])
        for idx in pc_indices
    ]


def euclidean(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def plot_pc(samples, unmasked, masked, pc_index, sign, png_path, pdf_path):
    x = [unmasked[sid]["pcs"][pc_index] for sid in samples]
    y = [sign * masked[sid]["pcs"][pc_index] for sid in samples]

    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    ax.scatter(x, y, s=18, alpha=0.65, linewidths=0)

    lo = min(x + y)
    hi = max(x + y)
    pad = 0.04 * (hi - lo) if hi > lo else 1.0
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], linewidth=0.8)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel(f"Unmasked PC{pc_index + 1}")
    ax.set_ylabel(f"High-LD-masked PC{pc_index + 1} (sign-aligned)")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)


def plot_side_by_side_population_pca(
    samples,
    unmasked,
    masked,
    unmasked_eval,
    masked_eval,
    signs,
    x_index,
    y_index,
    png_path,
    pdf_path,
):
    # Align masked PC signs to the unmasked solution before comparing panels.
    unmasked_xy = [
        (
            unmasked[sid]["pcs"][x_index],
            unmasked[sid]["pcs"][y_index],
            unmasked[sid]["population"],
        )
        for sid in samples
    ]
    masked_xy = [
        (
            signs[x_index] * masked[sid]["pcs"][x_index],
            signs[y_index] * masked[sid]["pcs"][y_index],
            masked[sid]["population"],
        )
        for sid in samples
    ]

    all_x = [x for x, _, _ in unmasked_xy] + [x for x, _, _ in masked_xy]
    all_y = [y for _, y, _ in unmasked_xy] + [y for _, y, _ in masked_xy]
    x_lo, x_hi = min(all_x), max(all_x)
    y_lo, y_hi = min(all_y), max(all_y)
    x_pad = 0.04 * (x_hi - x_lo) if x_hi > x_lo else 1.0
    y_pad = 0.04 * (y_hi - y_lo) if y_hi > y_lo else 1.0

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.5), sharex=True, sharey=True)
    panels = [
        ("Unmasked", unmasked_xy, unmasked_eval),
        ("High-LD masked", masked_xy, masked_eval),
    ]

    for ax, (title, values, evals) in zip(axes, panels):
        for pop in POP_ORDER:
            subset = [(x, y) for x, y, p in values if p == pop]
            xs = [x for x, _ in subset]
            ys = [y for _, y in subset]
            if pop == "CT":
                ax.scatter(
                    xs,
                    ys,
                    label="Cinque Terre",
                    marker="*",
                    s=72,
                    linewidths=0.7,
                    edgecolors="black",
                    zorder=5,
                )
            else:
                ax.scatter(
                    xs,
                    ys,
                    label=pop,
                    marker="o",
                    s=24,
                    alpha=0.72,
                    linewidths=0,
                )

        ax.axhline(0, linewidth=0.5, alpha=0.35)
        ax.axvline(0, linewidth=0.5, alpha=0.35)
        ax.set_xlim(x_lo - x_pad, x_hi + x_pad)
        ax.set_ylim(y_lo - y_pad, y_hi + y_pad)
        ax.set_title(title)
        ax.set_xlabel(
            f"PC{x_index + 1} ({evals[x_index + 1]['variance_percent']:.2f}%)"
        )

    axes[0].set_ylabel(
        f"PC{y_index + 1} ({unmasked_eval[y_index + 1]['variance_percent']:.2f}%)"
    )
    axes[1].set_ylabel(
        f"PC{y_index + 1} ({masked_eval[y_index + 1]['variance_percent']:.2f}%)"
    )
    axes[0].legend(frameon=False, fontsize=9, ncol=2)

    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--unmasked-coordinates", required=True)
    parser.add_argument("--masked-coordinates", required=True)
    parser.add_argument("--unmasked-eigenvalues", required=True)
    parser.add_argument("--masked-eigenvalues", required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    parser.add_argument("--n-components", type=int, required=True)
    parser.add_argument("--correlations-out", required=True)
    parser.add_argument("--correlation-matrix-out", required=True)
    parser.add_argument("--centroid-shifts-out", required=True)
    parser.add_argument("--ct-reference-distances-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--pc1-png", required=True)
    parser.add_argument("--pc1-pdf", required=True)
    parser.add_argument("--pc2-png", required=True)
    parser.add_argument("--pc2-pdf", required=True)
    parser.add_argument("--pc3-png", required=True)
    parser.add_argument("--pc3-pdf", required=True)
    parser.add_argument("--side-by-side-pc12-png", required=True)
    parser.add_argument("--side-by-side-pc12-pdf", required=True)
    parser.add_argument("--side-by-side-pc23-png", required=True)
    parser.add_argument("--side-by-side-pc23-pdf", required=True)
    args = parser.parse_args()

    unmasked = read_coordinates(Path(args.unmasked_coordinates), args.n_components)
    masked = read_coordinates(Path(args.masked_coordinates), args.n_components)
    unmasked_eval = read_eigenvalues(Path(args.unmasked_eigenvalues), args.n_components)
    masked_eval = read_eigenvalues(Path(args.masked_eigenvalues), args.n_components)

    if set(unmasked) != set(masked):
        raise SystemExit("ERROR: masked and unmasked PCA sample sets differ")
    samples = sorted(unmasked)
    if len(samples) != args.expected_samples:
        raise SystemExit(
            f"ERROR: expected {args.expected_samples} matched samples, found {len(samples)}"
        )
    for sid in samples:
        if unmasked[sid]["population"] != masked[sid]["population"]:
            raise SystemExit(f"ERROR: population label differs between panels for {sid}")

    correlation_rows = []
    signs = []
    matrix = []
    for i in range(args.n_components):
        ux = [unmasked[sid]["pcs"][i] for sid in samples]
        row = []
        for j in range(args.n_components):
            my = [masked[sid]["pcs"][j] for sid in samples]
            row.append(pearson(ux, my))
        matrix.append(row)

        same_r = row[i]
        sign = 1.0 if same_r >= 0 else -1.0
        signs.append(sign)
        aligned = [sign * masked[sid]["pcs"][i] for sid in samples]
        same_abs = abs(same_r)
        aligned_rmse = rmse(ux, aligned)
        best_j = max(range(args.n_components), key=lambda j: abs(row[j]))
        correlation_rows.append(
            [
                i + 1,
                same_r,
                same_abs,
                int(sign),
                aligned_rmse,
                best_j + 1,
                row[best_j],
                abs(row[best_j]),
                unmasked_eval[i + 1]["eigenvalue"],
                masked_eval[i + 1]["eigenvalue"],
                unmasked_eval[i + 1]["variance_percent"],
                masked_eval[i + 1]["variance_percent"],
            ]
        )

    outdir = Path(args.summary_out).parent
    outdir.mkdir(parents=True, exist_ok=True)

    with Path(args.correlations_out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "PC",
                "pearson_r_same_pc",
                "abs_pearson_r_same_pc",
                "masked_sign_alignment",
                "aligned_rmse",
                "best_matching_masked_PC",
                "best_match_r",
                "best_match_abs_r",
                "unmasked_eigenvalue",
                "masked_eigenvalue",
                "unmasked_variance_percent",
                "masked_variance_percent",
            ]
        )
        writer.writerows(correlation_rows)

    with Path(args.correlation_matrix_out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["unmasked_PC"] + [f"masked_PC{i}" for i in range(1, args.n_components + 1)])
        for i, row in enumerate(matrix, start=1):
            writer.writerow([f"PC{i}"] + row)

    population_samples = {
        pop: [sid for sid in samples if unmasked[sid]["population"] == pop]
        for pop in POP_ORDER
    }
    if any(not population_samples[pop] for pop in POP_ORDER):
        missing = [pop for pop in POP_ORDER if not population_samples[pop]]
        raise SystemExit(f"ERROR: empty expected population(s): {missing}")

    pc_indices = [0, 1, 2]
    with Path(args.centroid_shifts_out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "population",
                "n_samples",
                "unmasked_PC1_mean",
                "unmasked_PC2_mean",
                "unmasked_PC3_mean",
                "masked_PC1_mean_aligned",
                "masked_PC2_mean_aligned",
                "masked_PC3_mean_aligned",
                "centroid_shift_PC1_PC3",
            ]
        )
        for pop in POP_ORDER:
            ids = population_samples[pop]
            u = centroid(unmasked, ids, pc_indices)
            m_raw = centroid(masked, ids, pc_indices)
            m = [m_raw[k] * signs[k] for k in range(3)]
            writer.writerow([pop, len(ids)] + u + m + [euclidean(u, m)])

    ct_ids = population_samples["CT"]
    ct_u = centroid(unmasked, ct_ids, pc_indices)
    ct_m_raw = centroid(masked, ct_ids, pc_indices)
    ct_m = [ct_m_raw[k] * signs[k] for k in range(3)]

    distance_rows = []
    for pop in REF_POPS:
        ids = population_samples[pop]
        ref_u = centroid(unmasked, ids, pc_indices)
        ref_m_raw = centroid(masked, ids, pc_indices)
        ref_m = [ref_m_raw[k] * signs[k] for k in range(3)]
        du = euclidean(ct_u, ref_u)
        dm = euclidean(ct_m, ref_m)
        distance_rows.append([pop, du, dm, dm - du])

    with Path(args.ct_reference_distances_out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "reference_population",
                "ct_centroid_distance_unmasked_PC1_PC3",
                "ct_centroid_distance_masked_PC1_PC3",
                "masked_minus_unmasked_distance",
            ]
        )
        writer.writerows(distance_rows)

    abs_same = [row[2] for row in correlation_rows]
    first3 = abs_same[:3]
    same_pc_best_matches = sum(
        1 for i, row in enumerate(correlation_rows, start=1) if row[5] == i
    )
    max_centroid_shift = 0.0
    for pop in POP_ORDER:
        ids = population_samples[pop]
        u = centroid(unmasked, ids, pc_indices)
        m_raw = centroid(masked, ids, pc_indices)
        m = [m_raw[k] * signs[k] for k in range(3)]
        max_centroid_shift = max(max_centroid_shift, euclidean(u, m))

    with Path(args.summary_out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["matched_samples", len(samples)])
        writer.writerow(["components_compared", args.n_components])
        writer.writerow(["mean_abs_same_pc_correlation_first3", mean(first3)])
        writer.writerow(["minimum_abs_same_pc_correlation_first3", min(first3)])
        writer.writerow(["mean_abs_same_pc_correlation_first10", mean(abs_same)])
        writer.writerow(["minimum_abs_same_pc_correlation_first10", min(abs_same)])
        writer.writerow(["same_index_is_best_match_count_first10", same_pc_best_matches])
        writer.writerow(["max_population_centroid_shift_PC1_PC3", max_centroid_shift])
        for i in range(3):
            writer.writerow([f"PC{i+1}_abs_correlation", abs_same[i]])
            writer.writerow([f"PC{i+1}_aligned_rmse", correlation_rows[i][4]])

    plot_pc(samples, unmasked, masked, 0, signs[0], args.pc1_png, args.pc1_pdf)
    plot_pc(samples, unmasked, masked, 1, signs[1], args.pc2_png, args.pc2_pdf)
    plot_pc(samples, unmasked, masked, 2, signs[2], args.pc3_png, args.pc3_pdf)

    plot_side_by_side_population_pca(
        samples,
        unmasked,
        masked,
        unmasked_eval,
        masked_eval,
        signs,
        0,
        1,
        args.side_by_side_pc12_png,
        args.side_by_side_pc12_pdf,
    )
    plot_side_by_side_population_pca(
        samples,
        unmasked,
        masked,
        unmasked_eval,
        masked_eval,
        signs,
        1,
        2,
        args.side_by_side_pc23_png,
        args.side_by_side_pc23_pdf,
    )


if __name__ == "__main__":
    main()
