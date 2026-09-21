# ============================================================
# Individual pre-QC: diagnostic autosomal sample QC
# ============================================================
# This stage screens the 50 delivered samples before any exclusion.
# It is intentionally diagnostic: flags identify samples for review and do
# NOT remove samples automatically.
#
# Primary technical metrics:
#   - average depth
#   - genotype missingness / call rate
#   - heterozygosity
# Secondary diagnostics:
#   - non-reference variant count
#   - singleton fraction
#   - per-sample Ti/Tv
#
# Statistics are restricted to autosomes 1-22 so sex chromosomes do not
# introduce sex-dependent differences into the sample comparison.

AUTOSOMES = ",".join(str(chrom) for chrom in range(1, 23))


rule sample_preqc_bcftools_stats:
    input:
        vcf=lambda wildcards: config["local_paths"]["target_vcf"],
        preqc="results/preqc/preqc_done.txt"
    output:
        stats="results/sample_preqc/autosomal.per_sample.bcftools.stats.txt"
    log:
        "logs/sample_preqc/bcftools_stats.log"
    params:
        regions=AUTOSOMES
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/sample_preqc logs/sample_preqc
        bcftools stats -s - -r {params.regions} {input.vcf} \
            > {output.stats} 2> {log}
        """


rule sample_preqc_review_flags:
    input:
        stats="results/sample_preqc/autosomal.per_sample.bcftools.stats.txt"
    output:
        metrics="results/sample_preqc/sample_metrics_autosomal.tsv",
        flags="results/sample_preqc/sample_review_flags.tsv",
        summary="results/sample_preqc/sample_preqc_summary.tsv"
    params:
        robust_z_threshold=config["sample_preqc"]["robust_z_threshold"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python workflow/exploratory/sample_preqc_autosomal.py \
            --stats {input.stats} \
            --metrics-out {output.metrics} \
            --flags-out {output.flags} \
            --summary-out {output.summary} \
            --robust-z-threshold {params.robust_z_threshold}
        """
