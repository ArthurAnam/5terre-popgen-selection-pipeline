#!/usr/bin/env python3
import csv
from pathlib import Path
import matplotlib.pyplot as plt

EXCLUDED={"TSBC6060","TSBC6389","TSBD8047","TSBD8199"}

def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def main():
    metrics=read_tsv("results/sample_preqc/sample_metrics_autosomal.tsv")
    geno=read_tsv("results/sample_preqc/sample_genotype_field_diagnostics.tsv")
    geno_by={r["sample"]:r for r in geno}
    outdir=Path("results/sample_preqc/plots")
    outdir.mkdir(parents=True, exist_ok=True)

    fig,ax=plt.subplots(figsize=(7,5))
    for r in metrics:
        s=r["sample"]
        x=float(r["average_depth"])
        y=float(r["missing_rate"])*100
        ax.scatter(x,y,s=30)
        if s in EXCLUDED:
            ax.annotate(s,(x,y),xytext=(5,5),textcoords="offset points",fontsize=8)
    ax.set_xlabel("Average depth")
    ax.set_ylabel("Missing genotypes (%)")
    ax.set_title("Individual pre-QC: depth vs missingness")
    fig.tight_layout()
    fig.savefig(outdir/"depth_vs_missingness.png",dpi=180)
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(7,5))
    for r in metrics:
        s=r["sample"]
        if s not in geno_by:
            continue
        x=float(r["het_rate"])
        y=float(geno_by[s]["mean_het_abs_balance_deviation"])
        ax.scatter(x,y,s=30)
        if s in EXCLUDED:
            ax.annotate(s,(x,y),xytext=(5,5),textcoords="offset points",fontsize=8)
    ax.set_xlabel("Heterozygosity rate")
    ax.set_ylabel("Mean |heterozygous allele balance - 0.5|")
    ax.set_title("Individual pre-QC: heterozygosity vs allelic imbalance")
    fig.tight_layout()
    fig.savefig(outdir/"heterozygosity_vs_allelic_imbalance.png",dpi=180)
    plt.close(fig)

if __name__=="__main__":
    main()
