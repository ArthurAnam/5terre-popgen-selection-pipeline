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

    readme = """Individual sample pre-QC plots
================================

General
-------
Each point represents one of the 50 WGS samples in the delivered post-VQSR VCF.
All sample-level metrics used here are restricted to autosomes 1-22.

Only the four samples selected for exclusion after manual QC review are labelled
with their sample IDs:
  TSBC6060
  TSBC6389
  TSBD8047
  TSBD8199

Unlabelled points are the other samples. Point colours are purely graphical and
do NOT encode population, QC status, sex, batch, or any other biological group.

1) depth_vs_missingness.png
---------------------------
X axis: Average depth
  Per-sample average sequencing depth reported by bcftools stats on autosomes
  1-22.

Y axis: Missing genotypes (%)
  Percentage of autosomal sites with a missing genotype for that sample:
      100 * n_missing / (n_called + n_missing)

Interpretation:
  Samples far to the left have lower average depth.
  Samples higher on the plot have more missing genotype calls.
  TSBD8199 is visually separated from the main cluster because it combines
  very low depth with markedly increased missingness.

2) heterozygosity_vs_allelic_imbalance.png
------------------------------------------
X axis: Heterozygosity rate
  Fraction of called autosomal genotypes that are heterozygous:
      n_heterozygous / n_called

Y axis: Mean |heterozygous allele balance - 0.5|
  For heterozygous genotypes with AD available, allele balance is calculated as:
      ALT depth / (REF depth + ALT depth)
  The plotted value is the sample mean of the absolute deviation of this
  balance from 0.5.

Interpretation:
  Larger X values indicate more heterozygous genotype calls.
  Larger Y values indicate stronger average allelic imbalance among
  heterozygous calls.
  TSBC6060, TSBC6389 and TSBD8047 combine genome-wide excess
  heterozygosity with abnormal allelic balance relative to the cohort.

Important
---------
These plots are diagnostic summaries. They are not stand-alone automatic
filtering rules. Final sample exclusions were based on the combined evidence
from autosomal sample metrics, chromosome-by-chromosome consistency, and
GT/DP/GQ/AD diagnostics.
"""
    (outdir/"README.txt").write_text(readme, encoding="utf-8")

if __name__=="__main__":
    main()
