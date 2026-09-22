#!/usr/bin/env python3
import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

def mean(x): return statistics.mean(x) if x else 0.0
def median(x): return statistics.median(x) if x else 0.0

def pearson(xs, ys):
    if len(xs) < 2: return float("nan")
    mx, my = mean(xs), mean(ys)
    dx=[x-mx for x in xs]; dy=[y-my for y in ys]
    den=math.sqrt(sum(v*v for v in dx)*sum(v*v for v in dy))
    return sum(x*y for x,y in zip(dx,dy))/den if den else float("nan")

def ranks(vals):
    order=sorted(enumerate(vals), key=lambda z:z[1]); out=[0.0]*len(vals); i=0
    while i<len(order):
        j=i+1
        while j<len(order) and order[j][1]==order[i][1]: j+=1
        r=(i+1+j)/2.0
        for k in range(i,j): out[order[k][0]]=r
        i=j
    return out

def spearman(xs,ys): return pearson(ranks(xs),ranks(ys))

def read_annotations(path):
    out={}
    with open(path,encoding="utf-8") as h:
        r=csv.DictReader(h,delimiter="\t")
        for row in r: out[row["sample_id"]]=row["population"]
    return out

def read_fam(path):
    rows=[]
    with open(path,encoding="utf-8") as h:
        for line in h:
            if line.strip():
                f=line.split(); rows.append((f[0],f[1]))
    return rows

def read_hom(path):
    d=defaultdict(list); header=None
    with open(path,encoding="utf-8") as h:
        for line in h:
            if not line.strip(): continue
            f=line.split()
            if header is None:
                header=f
                if not {"IID","KB"}.issubset(header): raise SystemExit(f"ERROR: bad .hom header: {path}")
                continue
            row=dict(zip(header,f)); d[row["IID"]].append(float(row["KB"]))
    return d

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--joint-hom",required=True); p.add_argument("--allpop-hom",required=True)
    p.add_argument("--fam",required=True); p.add_argument("--annotations",required=True)
    p.add_argument("--expected-samples",type=int,required=True)
    p.add_argument("--individual-out",required=True); p.add_argument("--population-out",required=True)
    p.add_argument("--comparison-out",required=True); p.add_argument("--correlations-out",required=True)
    p.add_argument("--readme-out",required=True)
    a=p.parse_args()

    ann=read_annotations(a.annotations); samples=read_fam(a.fam)
    if len(samples)!=a.expected_samples: raise SystemExit("ERROR: unexpected sample count")
    data={"joint_maf0.05":read_hom(a.joint_hom),
          "all_population_maf0.05_intersection":read_hom(a.allpop_hom)}

    rows=[]; lookup={}
    for cond,hom in data.items():
        for fid,iid in samples:
            lengths=hom.get(iid,[])
            rec={"condition":cond,"population":ann[iid],"FID":fid,"IID":iid,
                 "N_ROH":len(lengths),"TOTAL_ROH_MB":sum(lengths)/1000.0,
                 "MEAN_ROH_MB":mean(lengths)/1000.0,
                 "MEDIAN_ROH_MB":median(lengths)/1000.0,
                 "MAX_ROH_MB":max(lengths,default=0.0)/1000.0}
            rows.append(rec); lookup[(cond,iid)]=rec

    Path(a.individual_out).parent.mkdir(parents=True,exist_ok=True)
    fields=["condition","population","FID","IID","N_ROH","TOTAL_ROH_MB","MEAN_ROH_MB","MEDIAN_ROH_MB","MAX_ROH_MB"]
    with open(a.individual_out,"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=fields); w.writeheader(); w.writerows(rows)

    grouped=defaultdict(list)
    for r in rows: grouped[(r["condition"],r["population"])].append(r)
    pfields=["condition","population","n_samples","mean_N_ROH","median_N_ROH","mean_TOTAL_ROH_MB","median_TOTAL_ROH_MB","mean_MAX_ROH_MB","median_MAX_ROH_MB"]
    with open(a.population_out,"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=pfields); w.writeheader()
        for cond in data:
            for pop in sorted(set(ann.values())):
                rs=grouped[(cond,pop)]
                w.writerow({"condition":cond,"population":pop,"n_samples":len(rs),
                    "mean_N_ROH":mean([x["N_ROH"] for x in rs]),"median_N_ROH":median([x["N_ROH"] for x in rs]),
                    "mean_TOTAL_ROH_MB":mean([x["TOTAL_ROH_MB"] for x in rs]),"median_TOTAL_ROH_MB":median([x["TOTAL_ROH_MB"] for x in rs]),
                    "mean_MAX_ROH_MB":mean([x["MAX_ROH_MB"] for x in rs]),"median_MAX_ROH_MB":median([x["MAX_ROH_MB"] for x in rs])})

    cfields=["population","FID","IID","joint_N_ROH","allpop_N_ROH","joint_TOTAL_ROH_MB","allpop_TOTAL_ROH_MB",
             "delta_TOTAL_ROH_MB_allpop_minus_joint","retained_TOTAL_ROH_fraction","joint_MAX_ROH_MB","allpop_MAX_ROH_MB"]
    with open(a.comparison_out,"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=cfields); w.writeheader()
        for fid,iid in samples:
            x=lookup[("joint_maf0.05",iid)]; y=lookup[("all_population_maf0.05_intersection",iid)]
            frac=y["TOTAL_ROH_MB"]/x["TOTAL_ROH_MB"] if x["TOTAL_ROH_MB"] else (1.0 if y["TOTAL_ROH_MB"]==0 else float("nan"))
            w.writerow({"population":ann[iid],"FID":fid,"IID":iid,"joint_N_ROH":x["N_ROH"],"allpop_N_ROH":y["N_ROH"],
                "joint_TOTAL_ROH_MB":x["TOTAL_ROH_MB"],"allpop_TOTAL_ROH_MB":y["TOTAL_ROH_MB"],
                "delta_TOTAL_ROH_MB_allpop_minus_joint":y["TOTAL_ROH_MB"]-x["TOTAL_ROH_MB"],
                "retained_TOTAL_ROH_fraction":frac,"joint_MAX_ROH_MB":x["MAX_ROH_MB"],"allpop_MAX_ROH_MB":y["MAX_ROH_MB"]})

    corrfields=["population","metric","n_samples","pearson","spearman"]
    with open(a.correlations_out,"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=corrfields); w.writeheader()
        for pop in ["ALL"]+sorted(set(ann.values())):
            ids=[iid for _,iid in samples if pop=="ALL" or ann[iid]==pop]
            for metric in ["N_ROH","TOTAL_ROH_MB","MAX_ROH_MB"]:
                xs=[lookup[("joint_maf0.05",iid)][metric] for iid in ids]
                ys=[lookup[("all_population_maf0.05_intersection",iid)][metric] for iid in ids]
                w.writerow({"population":pop,"metric":metric,"n_samples":len(ids),
                            "pearson":pearson(xs,ys),"spearman":spearman(xs,ys)})

    with open(a.readme_out,"w",encoding="utf-8") as h:
        h.write("ROH MAF-scope call sensitivity\n==============================\n\n")
        h.write("Both calls use identical ROH parameters: >=1.5 Mb, 50 SNP/window, 50 SNP/run, ")
        h.write("density 50 kb/SNP, gap 500 kb, 5 missing/window, 1 heterozygote/window, threshold 0.05.\n")
        h.write("Only marker ascertainment differs: joint MAF>=0.05 versus MAF>=0.05 within every population.\n")

if __name__=="__main__":
    main()
