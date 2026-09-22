#!/usr/bin/env python3
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

def read_annotations(path):
    out = {}
    with open(path, "r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        if not {"sample_id", "population"}.issubset(r.fieldnames or []):
            raise SystemExit("ERROR: annotations need sample_id and population")
        for row in r:
            if row["sample_id"] in out:
                raise SystemExit(f"ERROR: duplicate annotation {row['sample_id']}")
            out[row["sample_id"]] = row["population"]
    return out

def read_fam(path):
    rows=[]; seen=set()
    with open(path, "r", encoding="utf-8") as h:
        for line in h:
            if not line.strip(): continue
            f=line.split()
            if len(f) < 2: raise SystemExit("ERROR: malformed FAM")
            fid,iid=f[0],f[1]
            if iid in seen: raise SystemExit(f"ERROR: duplicate IID {iid}")
            seen.add(iid); rows.append((fid,iid))
    return rows

def read_hom(path, min_kb):
    rows=[]
    with open(path, "r", encoding="utf-8") as h:
        header=None
        for line in h:
            if not line.strip(): continue
            f=line.split()
            if header is None:
                header=f
                required={"FID","IID","CHR","SNP1","SNP2","POS1","POS2","KB","NSNP","DENSITY","PHOM","PHET"}
                if not required.issubset(set(header)):
                    raise SystemExit(f"ERROR: unexpected .hom header in {path}: {header}")
                continue
            row=dict(zip(header,f))
            kb=float(row["KB"])
            if kb < min_kb: continue
            rows.append({
                "FID":row["FID"],"IID":row["IID"],"CHR":row["CHR"],
                "SNP1":row["SNP1"],"SNP2":row["SNP2"],
                "POS1":int(row["POS1"]),"POS2":int(row["POS2"]),
                "KB":kb,"NSNP":int(row["NSNP"]),
                "DENSITY":float(row["DENSITY"]),
                "PHOM":float(row["PHOM"]),"PHET":float(row["PHET"])
            })
    return rows

def mean(x): return statistics.mean(x) if x else 0.0
def median(x): return statistics.median(x) if x else 0.0

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--primary-hom",required=True)
    p.add_argument("--vif-hom",required=True)
    p.add_argument("--fam",required=True)
    p.add_argument("--annotations",required=True)
    p.add_argument("--min-kb",type=float,required=True)
    p.add_argument("--expected-samples",type=int,required=True)
    p.add_argument("--segments-out",required=True)
    p.add_argument("--individual-out",required=True)
    p.add_argument("--population-out",required=True)
    p.add_argument("--comparison-out",required=True)
    p.add_argument("--readme-out",required=True)
    a=p.parse_args()

    ann=read_annotations(a.annotations)
    samples=read_fam(a.fam)
    if len(samples) != a.expected_samples:
        raise SystemExit(f"ERROR: expected {a.expected_samples} samples, found {len(samples)}")
    missing=[iid for _,iid in samples if iid not in ann]
    if missing: raise SystemExit(f"ERROR: {len(missing)} samples lack annotations")

    conditions={
        "population_history_unpruned": read_hom(a.primary_hom,a.min_kb),
        "light_vif_howrigan": read_hom(a.vif_hom,a.min_kb),
    }
    by_sample={}
    for cond,rows in conditions.items():
        d=defaultdict(list)
        for r in rows:
            if r["IID"] not in ann:
                raise SystemExit(f"ERROR: ROH sample {r['IID']} lacks annotation")
            d[r["IID"]].append(r)
        by_sample[cond]=d

    Path(a.segments_out).parent.mkdir(parents=True,exist_ok=True)
    with open(a.segments_out,"w",encoding="utf-8",newline="") as h:
        fields=["condition","population","FID","IID","CHR","SNP1","SNP2","POS1","POS2","KB","NSNP","DENSITY_KB_PER_SNP","PHOM","PHET"]
        w=csv.DictWriter(h,delimiter="\t",fieldnames=fields); w.writeheader()
        for cond,rows in conditions.items():
            for r in rows:
                w.writerow({
                    "condition":cond,"population":ann[r["IID"]],
                    "FID":r["FID"],"IID":r["IID"],"CHR":r["CHR"],
                    "SNP1":r["SNP1"],"SNP2":r["SNP2"],"POS1":r["POS1"],"POS2":r["POS2"],
                    "KB":r["KB"],"NSNP":r["NSNP"],"DENSITY_KB_PER_SNP":r["DENSITY"],
                    "PHOM":r["PHOM"],"PHET":r["PHET"]
                })

    indiv=[]
    for cond in conditions:
        for fid,iid in samples:
            lengths=[r["KB"] for r in by_sample[cond].get(iid,[])]
            indiv.append({
                "condition":cond,"population":ann[iid],"FID":fid,"IID":iid,
                "N_ROH":len(lengths),
                "TOTAL_ROH_MB":sum(lengths)/1000.0,
                "MEAN_ROH_MB":mean(lengths)/1000.0,
                "MEDIAN_ROH_MB":median(lengths)/1000.0,
                "MAX_ROH_MB":max(lengths,default=0.0)/1000.0,
            })
    indiv_fields=["condition","population","FID","IID","N_ROH","TOTAL_ROH_MB","MEAN_ROH_MB","MEDIAN_ROH_MB","MAX_ROH_MB"]
    with open(a.individual_out,"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=indiv_fields); w.writeheader(); w.writerows(indiv)

    grouped=defaultdict(list)
    for r in indiv: grouped[(r["condition"],r["population"])].append(r)
    pop_fields=["condition","population","n_samples","mean_N_ROH","median_N_ROH","mean_TOTAL_ROH_MB","median_TOTAL_ROH_MB","mean_MAX_ROH_MB","median_MAX_ROH_MB"]
    with open(a.population_out,"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=pop_fields); w.writeheader()
        for cond in conditions:
            for pop in sorted(set(ann.values())):
                rows=grouped[(cond,pop)]
                n=[r["N_ROH"] for r in rows]; total=[r["TOTAL_ROH_MB"] for r in rows]; mx=[r["MAX_ROH_MB"] for r in rows]
                w.writerow({
                    "condition":cond,"population":pop,"n_samples":len(rows),
                    "mean_N_ROH":mean(n),"median_N_ROH":median(n),
                    "mean_TOTAL_ROH_MB":mean(total),"median_TOTAL_ROH_MB":median(total),
                    "mean_MAX_ROH_MB":mean(mx),"median_MAX_ROH_MB":median(mx)
                })

    lookup={(r["condition"],r["IID"]):r for r in indiv}
    cmp_fields=["population","FID","IID","primary_N_ROH","light_vif_N_ROH","delta_N_ROH_light_vif_minus_primary","primary_TOTAL_ROH_MB","light_vif_TOTAL_ROH_MB","delta_TOTAL_ROH_MB_light_vif_minus_primary","primary_MAX_ROH_MB","light_vif_MAX_ROH_MB"]
    with open(a.comparison_out,"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=cmp_fields); w.writeheader()
        for fid,iid in samples:
            x=lookup[("population_history_unpruned",iid)]
            y=lookup[("light_vif_howrigan",iid)]
            w.writerow({
                "population":ann[iid],"FID":fid,"IID":iid,
                "primary_N_ROH":x["N_ROH"],"light_vif_N_ROH":y["N_ROH"],
                "delta_N_ROH_light_vif_minus_primary":y["N_ROH"]-x["N_ROH"],
                "primary_TOTAL_ROH_MB":x["TOTAL_ROH_MB"],"light_vif_TOTAL_ROH_MB":y["TOTAL_ROH_MB"],
                "delta_TOTAL_ROH_MB_light_vif_minus_primary":y["TOTAL_ROH_MB"]-x["TOTAL_ROH_MB"],
                "primary_MAX_ROH_MB":x["MAX_ROH_MB"],"light_vif_MAX_ROH_MB":y["MAX_ROH_MB"]
            })

    with open(a.readme_out,"w",encoding="utf-8") as h:
        h.write(
            "ROH call sensitivity comparison\n"
            "===============================\n\n"
            f"Only ROH >= {a.min_kb/1000.0:g} Mb are compared.\n"
            "Primary: unpruned population-history framework.\n"
            "Sensitivity: light VIF pruning plus Howrigan-derived PLINK settings.\n"
            "FROH is deliberately deferred until the ROH definition and denominator are frozen.\n"
        )

if __name__ == "__main__":
    main()
