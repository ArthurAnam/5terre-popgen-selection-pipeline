#!/usr/bin/env python3
import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

GAPS = (100, 500, 1000)
HETS = (0, 1)
BASELINE = (1000, 1)

def read_annotations(path):
    mapping = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"sample_id", "population"}.issubset(reader.fieldnames or []):
            raise SystemExit("ERROR: annotations require sample_id and population")
        for row in reader:
            iid = row["sample_id"]
            if iid in mapping:
                raise SystemExit(f"ERROR: duplicate annotation for {iid}")
            mapping[iid] = row["population"]
    return mapping

def read_fam(path):
    rows = []
    seen = set()
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.split()
            if len(fields) < 2:
                raise SystemExit("ERROR: malformed FAM")
            fid, iid = fields[0], fields[1]
            if iid in seen:
                raise SystemExit(f"ERROR: duplicate IID in FAM: {iid}")
            seen.add(iid)
            rows.append((fid, iid))
    return rows

def read_hom(path):
    by_iid = defaultdict(list)
    with open(path, "r", encoding="utf-8") as handle:
        header = None
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.split()
            if header is None:
                header = fields
                required = {"FID", "IID", "KB"}
                if not required.issubset(set(header)):
                    raise SystemExit(f"ERROR: unexpected .hom header in {path}")
                continue
            row = dict(zip(header, fields))
            by_iid[row["IID"]].append(float(row["KB"]))
    return by_iid

def mean(values):
    return statistics.mean(values) if values else 0.0

def median(values):
    return statistics.median(values) if values else 0.0

def pearson(xs, ys):
    if len(xs) < 2:
        return float("nan")
    mx, my = mean(xs), mean(ys)
    dx = [x - mx for x in xs]
    dy = [y - my for y in ys]
    den = math.sqrt(sum(x*x for x in dx) * sum(y*y for y in dy))
    return sum(x*y for x, y in zip(dx, dy)) / den if den else float("nan")

def ranks(values):
    ordered = sorted(enumerate(values), key=lambda x: x[1])
    out = [0.0] * len(values)
    i = 0
    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and ordered[j][1] == ordered[i][1]:
            j += 1
        rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            out[ordered[k][0]] = rank
        i = j
    return out

def spearman(xs, ys):
    return pearson(ranks(xs), ranks(ys))

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True)
    p.add_argument("--fam", required=True)
    p.add_argument("--annotations", required=True)
    p.add_argument("--expected-samples", type=int, required=True)
    p.add_argument("--individual-out", required=True)
    p.add_argument("--population-out", required=True)
    p.add_argument("--comparison-out", required=True)
    p.add_argument("--correlations-out", required=True)
    p.add_argument("--readme-out", required=True)
    a = p.parse_args()

    samples = read_fam(a.fam)
    if len(samples) != a.expected_samples:
        raise SystemExit(f"ERROR: expected {a.expected_samples} samples, found {len(samples)}")
    ann = read_annotations(a.annotations)
    missing = [iid for _, iid in samples if iid not in ann]
    if missing:
        raise SystemExit(f"ERROR: {len(missing)} FAM samples lack annotations")

    summaries = {}
    rows = []
    for gap in GAPS:
        for het in HETS:
            label = f"gap{gap}_het{het}"
            hom = read_hom(Path(a.input_dir) / f"{label}.hom")
            lookup = {}
            for fid, iid in samples:
                lengths = hom.get(iid, [])
                rec = {
                    "condition": label,
                    "gap_kb": gap,
                    "window_het": het,
                    "population": ann[iid],
                    "FID": fid,
                    "IID": iid,
                    "N_ROH": len(lengths),
                    "TOTAL_ROH_MB": sum(lengths) / 1000.0,
                    "MEAN_ROH_MB": mean(lengths) / 1000.0,
                    "MEDIAN_ROH_MB": median(lengths) / 1000.0,
                    "MAX_ROH_MB": max(lengths, default=0.0) / 1000.0,
                }
                rows.append(rec)
                lookup[iid] = rec
            summaries[(gap, het)] = lookup

    indiv_fields = ["condition","gap_kb","window_het","population","FID","IID","N_ROH",
                    "TOTAL_ROH_MB","MEAN_ROH_MB","MEDIAN_ROH_MB","MAX_ROH_MB"]
    Path(a.individual_out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.individual_out, "w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, delimiter="\t", fieldnames=indiv_fields)
        w.writeheader()
        w.writerows(rows)

    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["condition"], row["population"])].append(row)

    pop_fields = ["condition","gap_kb","window_het","population","n_samples",
                  "mean_N_ROH","median_N_ROH","mean_TOTAL_ROH_MB","median_TOTAL_ROH_MB",
                  "mean_MAX_ROH_MB","median_MAX_ROH_MB"]
    with open(a.population_out, "w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, delimiter="\t", fieldnames=pop_fields)
        w.writeheader()
        for gap in GAPS:
            for het in HETS:
                label = f"gap{gap}_het{het}"
                for pop in sorted(set(ann.values())):
                    rs = grouped[(label, pop)]
                    n = [r["N_ROH"] for r in rs]
                    total = [r["TOTAL_ROH_MB"] for r in rs]
                    mx = [r["MAX_ROH_MB"] for r in rs]
                    w.writerow({
                        "condition": label, "gap_kb": gap, "window_het": het,
                        "population": pop, "n_samples": len(rs),
                        "mean_N_ROH": mean(n), "median_N_ROH": median(n),
                        "mean_TOTAL_ROH_MB": mean(total), "median_TOTAL_ROH_MB": median(total),
                        "mean_MAX_ROH_MB": mean(mx), "median_MAX_ROH_MB": median(mx),
                    })

    base = summaries[BASELINE]
    cmp_fields = ["condition","population","FID","IID","N_ROH","baseline_N_ROH",
                  "delta_N_ROH","TOTAL_ROH_MB","baseline_TOTAL_ROH_MB",
                  "delta_TOTAL_ROH_MB","MAX_ROH_MB","baseline_MAX_ROH_MB",
                  "delta_MAX_ROH_MB"]
    with open(a.comparison_out, "w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, delimiter="\t", fieldnames=cmp_fields)
        w.writeheader()
        for gap in GAPS:
            for het in HETS:
                label = f"gap{gap}_het{het}"
                cur = summaries[(gap, het)]
                for fid, iid in samples:
                    x, b = cur[iid], base[iid]
                    w.writerow({
                        "condition": label, "population": ann[iid], "FID": fid, "IID": iid,
                        "N_ROH": x["N_ROH"], "baseline_N_ROH": b["N_ROH"],
                        "delta_N_ROH": x["N_ROH"] - b["N_ROH"],
                        "TOTAL_ROH_MB": x["TOTAL_ROH_MB"],
                        "baseline_TOTAL_ROH_MB": b["TOTAL_ROH_MB"],
                        "delta_TOTAL_ROH_MB": x["TOTAL_ROH_MB"] - b["TOTAL_ROH_MB"],
                        "MAX_ROH_MB": x["MAX_ROH_MB"],
                        "baseline_MAX_ROH_MB": b["MAX_ROH_MB"],
                        "delta_MAX_ROH_MB": x["MAX_ROH_MB"] - b["MAX_ROH_MB"],
                    })

    corr_fields = ["condition","population","metric","n_samples","pearson","spearman"]
    with open(a.correlations_out, "w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, delimiter="\t", fieldnames=corr_fields)
        w.writeheader()
        populations = ["ALL"] + sorted(set(ann.values()))
        metric_map = {
            "N_ROH": "N_ROH",
            "TOTAL_ROH_MB": "TOTAL_ROH_MB",
            "MAX_ROH_MB": "MAX_ROH_MB",
        }
        for gap in GAPS:
            for het in HETS:
                label = f"gap{gap}_het{het}"
                cur = summaries[(gap, het)]
                for pop in populations:
                    ids = [iid for _, iid in samples if pop == "ALL" or ann[iid] == pop]
                    for metric, key in metric_map.items():
                        xs = [base[iid][key] for iid in ids]
                        ys = [cur[iid][key] for iid in ids]
                        w.writerow({
                            "condition": label, "population": pop, "metric": metric,
                            "n_samples": len(ids), "pearson": pearson(xs, ys),
                            "spearman": spearman(xs, ys),
                        })

    with open(a.readme_out, "w", encoding="utf-8") as handle:
        handle.write(
            "ROH gap x heterozygote sensitivity\n"
            "===================================\n\n"
            "Unpruned common MAF>=0.05 panel; ROH >=1.5 Mb; 50 SNP/window;\n"
            "50 SNP/run; density 50 kb/SNP; 5 missing/window; threshold 0.05.\n"
            "Only homozyg-gap (100, 500, 1000 kb) and homozyg-window-het (0, 1)\n"
            "are varied. gap1000_het1 is the current population-history baseline.\n"
            "This isolates the effects of gap bridging and heterozygote tolerance\n"
            "before the production ROH definition is frozen.\n"
        )

if __name__ == "__main__":
    main()
