#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path


def read_coords(path):
    out = {}
    with Path(path).open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        required = {"sample_id", "population", "PC1", "PC2"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"Missing columns in {path}: {sorted(missing)}")

        for row in r:
            sid = row["sample_id"]
            if sid in out:
                raise SystemExit(f"Duplicate sample: {sid}")
            out[sid] = (
                row["population"],
                float(row["PC1"]),
                float(row["PC2"]),
            )
    return out


def mean(x):
    return sum(x) / len(x)


def centre(points):
    mx = mean([x for x, _ in points])
    my = mean([y for _, y in points])
    return [(x - mx, y - my) for x, y in points], (mx, my)


def ss(points):
    return sum(x*x + y*y for x, y in points)


def optimal_orthogonal(x, y):
    # M = Y'X
    m00 = sum(yy[0] * xx[0] for xx, yy in zip(x, y))
    m01 = sum(yy[0] * xx[1] for xx, yy in zip(x, y))
    m10 = sum(yy[1] * xx[0] for xx, yy in zip(x, y))
    m11 = sum(yy[1] * xx[1] for xx, yy in zip(x, y))

    # det +1 rotation
    a = m00 + m11
    b = m10 - m01
    nr = math.hypot(a, b)
    if nr:
        cr, sr = a / nr, b / nr
    else:
        cr, sr = 1.0, 0.0
    score_rotation = nr

    # det -1 reflection
    a = m00 - m11
    b = m10 + m01
    nf = math.hypot(a, b)
    if nf:
        cf, sf = a / nf, b / nf
    else:
        cf, sf = 1.0, 0.0
    score_reflection = nf

    if score_rotation >= score_reflection:
        return (
            ((cr, -sr), (sr, cr)),
            1,
            score_rotation,
        )

    return (
        ((cf, sf), (sf, -cf)),
        -1,
        score_reflection,
    )


def transform(points, R):
    return [
        (
            x * R[0][0] + y * R[1][0],
            x * R[0][1] + y * R[1][1],
        )
        for x, y in points
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unmasked", required=True)
    ap.add_argument("--masked", required=True)
    ap.add_argument("--expected-samples", type=int, required=True)
    ap.add_argument("--summary-out", required=True)
    ap.add_argument("--displacements-out", required=True)
    args = ap.parse_args()

    u = read_coords(args.unmasked)
    m = read_coords(args.masked)

    if set(u) != set(m):
        raise SystemExit("Masked and unmasked sample sets differ")

    samples = sorted(u)

    if len(samples) != args.expected_samples:
        raise SystemExit(
            f"Expected {args.expected_samples} samples, found {len(samples)}"
        )

    for sid in samples:
        if u[sid][0] != m[sid][0]:
            raise SystemExit(f"Population mismatch for {sid}")

    X = [(u[s][1], u[s][2]) for s in samples]
    Y = [(m[s][1], m[s][2]) for s in samples]

    Xc, xmean = centre(X)
    Yc, ymean = centre(Y)

    sx = math.sqrt(ss(Xc))
    sy = math.sqrt(ss(Yc))
    if sx == 0 or sy == 0:
        raise SystemExit("Degenerate PCA configuration")

    # Symmetric standardized Procrustes disparity.
    Xn = [(x/sx, y/sx) for x, y in Xc]
    Yn = [(x/sy, y/sy) for x, y in Yc]

    Rn, detR, _ = optimal_orthogonal(Xn, Yn)
    YnR = transform(Yn, Rn)

    disparity = sum(
        (x1-x2)**2 + (y1-y2)**2
        for (x1, y1), (x2, y2) in zip(Xn, YnR)
    )

    normalized_rms_per_sample = math.sqrt(disparity / len(samples))

    # Full Procrustes fit in original unmasked PCA coordinate units:
    # translate masked to origin, rotate/reflect, then fit one isotropic scale.
    R, detR2, score = optimal_orthogonal(Xc, Yc)
    if detR2 != detR:
        raise SystemExit("Internal Procrustes orientation inconsistency")

    scale = score / ss(Yc)
    Yrot = transform(Yc, R)

    fitted = [
        (
            xmean[0] + scale*x,
            xmean[1] + scale*y,
        )
        for x, y in Yrot
    ]

    residuals = []
    rows = []

    for sid, xu, yf in zip(samples, X, fitted):
        dx = yf[0] - xu[0]
        dy = yf[1] - xu[1]
        d = math.hypot(dx, dy)
        residuals.append(d)
        rows.append([
            sid,
            u[sid][0],
            xu[0],
            xu[1],
            m[sid][1],
            m[sid][2],
            yf[0],
            yf[1],
            d,
        ])

    residuals_sorted = sorted(residuals)

    def quantile(q):
        pos = q * (len(residuals_sorted) - 1)
        lo = int(math.floor(pos))
        hi = int(math.ceil(pos))
        if lo == hi:
            return residuals_sorted[lo]
        f = pos - lo
        return (
            residuals_sorted[lo] * (1-f)
            + residuals_sorted[hi] * f
        )

    coord_rmse = math.sqrt(
        sum(d*d for d in residuals) / (2 * len(samples))
    )

    pop_values = {}
    pop_rows = {}

    for row in rows:
        pop = row[1]
        pop_values.setdefault(pop, []).append(row[-1])
        pop_rows.setdefault(pop, []).append(row)

    pop_centroid_shifts = {}

    for pop, prows in pop_rows.items():
        unmasked_centroid = (
            mean([r[2] for r in prows]),
            mean([r[3] for r in prows]),
        )
        aligned_masked_centroid = (
            mean([r[6] for r in prows]),
            mean([r[7] for r in prows]),
        )

        pop_centroid_shifts[pop] = math.hypot(
            aligned_masked_centroid[0] - unmasked_centroid[0],
            aligned_masked_centroid[1] - unmasked_centroid[1],
        )

    Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)

    with Path(args.displacements_out).open(
        "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "sample_id",
            "population",
            "unmasked_PC1",
            "unmasked_PC2",
            "masked_PC1",
            "masked_PC2",
            "masked_PC1_procrustes_aligned",
            "masked_PC2_procrustes_aligned",
            "euclidean_displacement",
        ])
        w.writerows(rows)

    with Path(args.summary_out).open(
        "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["matched_samples", len(samples)])
        w.writerow(["dimensions", 2])
        w.writerow(["standardized_procrustes_disparity", disparity])
        w.writerow([
            "normalized_rms_displacement_per_sample",
            normalized_rms_per_sample,
        ])
        w.writerow(["orthogonal_transform_determinant", detR])
        w.writerow(["fitted_isotropic_scale_masked_to_unmasked", scale])
        w.writerow(["coordinate_rmse_after_alignment", coord_rmse])
        w.writerow(["mean_sample_displacement", mean(residuals)])
        w.writerow(["median_sample_displacement", quantile(0.5)])
        w.writerow(["p95_sample_displacement", quantile(0.95)])
        w.writerow(["max_sample_displacement", max(residuals)])

        for pop in sorted(pop_values):
            vals = sorted(pop_values[pop])
            w.writerow([
                f"{pop}_mean_displacement",
                sum(vals) / len(vals),
            ])
            w.writerow([
                f"{pop}_max_displacement",
                max(vals),
            ])
            w.writerow([
                f"{pop}_centroid_shift_PC1_PC2",
                pop_centroid_shifts[pop],
            ])


if __name__ == "__main__":
    main()
