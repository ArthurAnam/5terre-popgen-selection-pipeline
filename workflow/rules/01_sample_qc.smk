# ============================================================
# Per-sample WGS QC inspired by Kars et al. supplementary QC
# ============================================================
# This stage is diagnostic: it identifies candidate sample outliers but does
# NOT remove samples. The four exclusion metrics mirror the paper:
# singletons, Ti/Tv, average depth, and total number of variants.
# Missing sites/rate are reported as additional diagnostics.

rule sample_qc_bcftools_stats:
    input:
        vcf=lambda wildcards: config["local_paths"]["target_vcf"],
        preqc="results/preqc/preqc_done.txt"
    output:
        stats="results/sample_qc/per_sample.bcftools.stats.txt"
    log:
        "logs/sample_qc/bcftools_stats.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/sample_qc logs/sample_qc
        bcftools stats -s - {input.vcf} > {output.stats} 2> {log}
        """


rule sample_qc_flag_outliers:
    input:
        stats="results/sample_qc/per_sample.bcftools.stats.txt"
    output:
        metrics="results/sample_qc/sample_metrics.tsv",
        outliers="results/sample_qc/sample_outliers_5MAD.tsv",
        summary="results/sample_qc/sample_qc_summary.tsv"
    params:
        mad_threshold=5.0
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python workflow/scripts/sample_qc_from_bcftools_stats.py \
            --stats {input.stats} \
            --metrics-out {output.metrics} \
            --outliers-out {output.outliers} \
            --summary-out {output.summary} \
            --mad-threshold {params.mad_threshold}
        """
