# ============================================================
# Genotype-level QC exploration and exact per-individual missingness
# ============================================================

GENOTYPE_QC_REGIONS_POST_SAMPLE_QC = ",".join(
    f"{chrom}:{start}-{start + 999999}"
    for chrom in range(1, 23)
    for start in (20000000, 35000000, 45000000)
)


rule explore_genotype_qc_and_individual_missingness:
    input:
        vcf="results/qc/00_post_sample_qc_input/cinque_terre.post_sample_qc.autosomes.vcf.gz",
        index="results/qc/00_post_sample_qc_input/cinque_terre.post_sample_qc.autosomes.vcf.gz.tbi",
        samples="results/qc/00_post_sample_qc_input/retained_samples_after_sample_qc.txt"
    output:
        genotype_summary="results/qc/02_genotype_qc_exploration/genotype_qc_metric_descriptive_statistics.tsv",
        genotype_plot="results/qc/02_genotype_qc_exploration/genotype_qc_metric_distributions.png",
        sample_missingness="results/qc/02_genotype_qc_exploration/individual_missingness_before_variant_filtering.tsv",
        readme="results/qc/02_genotype_qc_exploration/README.txt"
    log:
        "logs/qc/02_genotype_qc_exploration/genotype_qc_exploration.log"
    params:
        regions=GENOTYPE_QC_REGIONS_POST_SAMPLE_QC
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/02_genotype_qc_exploration logs/qc/02_genotype_qc_exploration

        n_samples=$(wc -l < {input.samples})

        {{
            printf 'CHROM\tPOS'
            while read -r sample; do
                printf '\t%s_GT\t%s_DP\t%s_GQ\t%s_AD' "$sample" "$sample" "$sample" "$sample"
            done < {input.samples}
            printf '\n'

            bcftools query                 -r {params.regions}                 -f '%CHROM\t%POS[\t%GT\t%DP\t%GQ\t%AD]\n'                 {input.vcf} 2> {log}
        }} | python workflow/scripts/explore_genotype_qc_metrics.py             --summary-out {output.genotype_summary}             --plot-out {output.genotype_plot}

        n_sites=$(bcftools index -n {input.vcf})

        bcftools stats -s - {input.vcf} 2>> {log}         | awk -v n_sites="$n_sites" -v cutoff="0.05" -F '\t' '
            BEGIN {
                OFS="\t";
                print "sample","n_sites","n_missing","missing_rate","above_0.05"
            }
            $1=="PSC" {
                rate=$14/n_sites;
                flag=(rate>cutoff ? "YES" : "NO");
                printf "%s\t%d\t%d\t%.8f\t%s\n", $3,n_sites,$14,rate,flag
            }
        ' > {output.sample_missingness}

        n_over=$(awk -F '\t' 'NR>1 && $5=="YES" {{n++}} END {{print n+0}}' {output.sample_missingness})

        cat > {output.readme} <<EOF
02_genotype_qc_exploration
==========================

Purpose
-------
Explore genotype-level quality after the four reviewed sample exclusions,
without applying DP, GQ or allele-balance hard filters.

Genotype-field exploration
--------------------------
FORMAT/DP, FORMAT/GQ and FORMAT/AD are summarized on the same deterministic
genome-wide subset used during sample pre-QC: three 1-Mb windows per autosome
(starting at 20, 35 and 45 Mb).

This subsampling is used only for descriptive distributions and keeps the
analysis reproducible and computationally manageable.

Individual missingness
----------------------
individual_missingness_before_variant_filtering.tsv is exact, not subsampled.
It is computed over all autosomal sites currently present in the 46-sample VCF.

The planned individual filter is >5% missingness. At this stage the table is
diagnostic only because the final individual-missingness decision will be
recomputed after the planned site-level variant filters.

Number of individuals currently above 5% missingness: $n_over

No genotype or individual is removed by this rule.
EOF
        """
