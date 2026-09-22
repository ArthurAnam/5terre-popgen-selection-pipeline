#!/usr/bin/env python3
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

MIN_KB = 1500.0
LONG_KB = 5000.0

def mean(values):
    return statistics.mean(values) if values else 0.0

def median(values):
    return statistics.median(values) if values else 0.0

def read_annotations(path):
    out = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"sample_id", "population"}.issubset(reader.fieldnames or []):
            raise SystemExit("ERROR: annotations require sample_id and population")
        for row in reader:
            iid = row["sample_id"]
            if iid in out:
                raise SystemExit(f"ERROR: duplicate annotation for {iid}")
            out[iid] = row["population"]
    return out

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
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        header = None
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.split()
            if header is None:
                header = fields
                required = {"FID","IID","CHR","SNP1","SNP2","POS1","POS2","KB","NSNP","DENSITY","PHOM","PHET"}
                if not required.issubset(set(header)):
                    raise SystemExit(f"ERROR: unexpected .hom header: {header}")
                continue
            row = dict(zip(header, fields))
            kb = float(row["KB"])
            if kb < MIN_KB:
                raise SystemExit(f"ERROR: production .hom contains ROH shorter than {MIN_KB} kb")
            rows.append({
                "FID": row["FID"], "IID": row["IID"], "CHR": row["CHR"],
                "SNP1": row["SNP1"], "SNP2": row["SNP2"],
                "POS1": int(row["POS1"]), "POS2": int(row["POS2"]),
                "KB": kb, "NSNP": int(row["NSNP"]), "DENSITY": float(row["DENSITY"]),
                "PHOM": float(row["PHOM"]), "PHET": float(row["PHET"])
            })
    return rows

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hom", required=True)
    p.add_argument("--fam", required=True)
    p.add_argument("--annotations", required=True)
    p.add_argument("--denominator-bp", type=int, required=True)
    p.add_argument("--expected-samples", type=int, required=True)
    p.add_argument("--segments-out", required=True)
    p.add_argument("--individual-out", required=True)
    p.add_argument("--population-out", required=True)
    p.add_argument("--readme-out", required=True)
    a = p.parse_args()

    ann = read_annotations(a.annotations)
    samples = read_fam(a.fam)
    if len(samples) != a.expected_samples:
        raise SystemExit(f"ERROR: expected {a.expected_samples} samples, found {len(samples)}")
    missing = [iid for _, iid in samples if iid not in ann]
    if missing:
        raise SystemExit(f"ERROR: {len(missing)} samples lack population annotation")

    segments = read_hom(a.hom)
    by_iid = defaultdict(list)
    for row in segments:
        if row["IID"] not in ann:
            raise SystemExit(f"ERROR: ROH IID {row['IID']} lacks annotation")
        by_iid[row["IID"]].append(row)

    Path(a.segments_out).parent.mkdir(parents=True, exist_ok=True)
    seg_fields = ["population","FID","IID","CHR","SNP1","SNP2","POS1","POS2","KB","MB","NSNP",
                  "DENSITY_KB_PER_SNP","PHOM","PHET","length_class"]
    with open(a.segments_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=seg_fields)
        writer.writeheader()
        for row in segments:
            writer.writerow({
                "population": ann[row["IID"]], "FID": row["FID"], "IID": row["IID"],
                "CHR": row["CHR"], "SNP1": row["SNP1"], "SNP2": row["SNP2"],
                "POS1": row["POS1"], "POS2": row["POS2"], "KB": row["KB"],
                "MB": row["KB"] / 1000.0, "NSNP": row["NSNP"],
                "DENSITY_KB_PER_SNP": row["DENSITY"], "PHOM": row["PHOM"], "PHET": row["PHET"],
                "length_class": "1.5_to_lt5_Mb" if row["KB"] < LONG_KB else "ge5_Mb"
            })

    individual = []
    for fid, iid in samples:
        segs = by_iid.get(iid, [])
        all_kb = [x["KB"] for x in segs]
        mid_kb = [x["KB"] for x in segs if x["KB"] < LONG_KB]
        long_kb = [x["KB"] for x in segs if x["KB"] >= LONG_KB]
        total_bp = sum(all_kb) * 1000.0
        mid_bp = sum(mid_kb) * 1000.0
        long_bp = sum(long_kb) * 1000.0
        individual.append({
            "population": ann[iid], "FID": fid, "IID": iid,
            "N_ROH_GE1_5MB": len(all_kb),
            "TOTAL_ROH_MB_GE1_5MB": sum(all_kb) / 1000.0,
            "FROH_GE1_5MB": total_bp / a.denominator_bp,
            "MEAN_ROH_MB_GE1_5MB": mean(all_kb) / 1000.0,
            "MEDIAN_ROH_MB_GE1_5MB": median(all_kb) / 1000.0,
            "MAX_ROH_MB": max(all_kb, default=0.0) / 1000.0,
            "N_ROH_1_5_TO_LT5MB": len(mid_kb),
            "TOTAL_ROH_MB_1_5_TO_LT5MB": sum(mid_kb) / 1000.0,
            "FROH_1_5_TO_LT5MB": mid_bp / a.denominator_bp,
            "N_ROH_GE5MB": len(long_kb),
            "TOTAL_ROH_MB_GE5MB": sum(long_kb) / 1000.0,
            "FROH_GE5MB": long_bp / a.denominator_bp,
        })

    indiv_fields = list(individual[0].keys())
    with open(a.individual_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=indiv_fields)
        writer.writeheader()
        writer.writerows(individual)

    grouped = defaultdict(list)
    for row in individual:
        grouped[row["population"]].append(row)

    pop_fields = ["population","n_samples",
                  "mean_N_ROH_GE1_5MB","median_N_ROH_GE1_5MB",
                  "mean_TOTAL_ROH_MB_GE1_5MB","median_TOTAL_ROH_MB_GE1_5MB",
                  "mean_FROH_GE1_5MB","median_FROH_GE1_5MB",
                  "mean_MAX_ROH_MB","median_MAX_ROH_MB",
                  "mean_N_ROH_GE5MB","median_N_ROH_GE5MB",
                  "mean_TOTAL_ROH_MB_GE5MB","median_TOTAL_ROH_MB_GE5MB",
                  "mean_FROH_GE5MB","median_FROH_GE5MB"]
    with open(a.population_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=pop_fields)
        writer.writeheader()
        for pop in sorted(grouped):
            rows = grouped[pop]
            writer.writerow({
                "population": pop, "n_samples": len(rows),
                "mean_N_ROH_GE1_5MB": mean([r["N_ROH_GE1_5MB"] for r in rows]),
                "median_N_ROH_GE1_5MB": median([r["N_ROH_GE1_5MB"] for r in rows]),
                "mean_TOTAL_ROH_MB_GE1_5MB": mean([r["TOTAL_ROH_MB_GE1_5MB"] for r in rows]),
                "median_TOTAL_ROH_MB_GE1_5MB": median([r["TOTAL_ROH_MB_GE1_5MB"] for r in rows]),
                "mean_FROH_GE1_5MB": mean([r["FROH_GE1_5MB"] for r in rows]),
                "median_FROH_GE1_5MB": median([r["FROH_GE1_5MB"] for r in rows]),
                "mean_MAX_ROH_MB": mean([r["MAX_ROH_MB"] for r in rows]),
                "median_MAX_ROH_MB": median([r["MAX_ROH_MB"] for r in rows]),
                "mean_N_ROH_GE5MB": mean([r["N_ROH_GE5MB"] for r in rows]),
                "median_N_ROH_GE5MB": median([r["N_ROH_GE5MB"] for r in rows]),
                "mean_TOTAL_ROH_MB_GE5MB": mean([r["TOTAL_ROH_MB_GE5MB"] for r in rows]),
                "median_TOTAL_ROH_MB_GE5MB": median([r["TOTAL_ROH_MB_GE5MB"] for r in rows]),
                "mean_FROH_GE5MB": mean([r["FROH_GE5MB"] for r in rows]),
                "median_FROH_GE5MB": median([r["FROH_GE5MB"] for r in rows]),
            })

    with open(a.readme_out, "w", encoding="utf-8") as handle:
        handle.write(
            "Production ROH and FROH\n"
            "=======================\n\n"
            "Primary marker panel: joint CT + 1000G EUR MAF >= 0.05, no LD pruning.\n"
            "PLINK ROH: 50-SNP scanning window; minimum 50 SNP per called ROH;\n"
            "minimum 1.5 Mb; density <=50 kb/SNP; gap <=500 kb; <=5 missing and\n"
            "<=1 heterozygous call per scanning window; hit threshold 0.05.\n"
            "No whole-run --homozyg-het cap is imposed.\n\n"
            f"FROH denominator: {a.denominator_bp} bp (2.77 Gb), the conventional\n"
            "human SNP-mappable autosomal distance.\n"
            "Primary FROH sums ROH >=1.5 Mb. Secondary length summaries split\n"
            "segments into 1.5-<5 Mb and >=5 Mb while retaining continuous lengths.\n"
        )

if __name__ == "__main__":
    main()
