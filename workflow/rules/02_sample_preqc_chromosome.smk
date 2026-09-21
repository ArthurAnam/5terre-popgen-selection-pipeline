# ============================================================
# Individual pre-QC: chromosome-level consistency
# ============================================================
# Checks whether sample-level anomalies are genome-wide or restricted to a
# small number of autosomes. This remains diagnostic and performs no filtering.

AUTOSOME_LIST = [str(chrom) for chrom in range(1, 23)]


rule sample_preqc_chr_bcftools_stats:
    input:
        vcf=lambda wildcards: config["local_paths"]["target_vcf"],
        preqc="results/preqc/preqc_done.txt"
    output:
        stats="results/sample_preqc/by_chromosome/chr{chrom}.bcftools.stats.txt"
    log:
        "logs/sample_preqc/by_chromosome/chr{chrom}.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/sample_preqc/by_chromosome logs/sample_preqc/by_chromosome

        bcftools stats -s - -r {wildcards.chrom} {input.vcf} \
            > {output.stats} 2> {log}
        """


rule sample_preqc_chr_consistency:
    input:
        stats=expand(
            "results/sample_preqc/by_chromosome/chr{chrom}.bcftools.stats.txt",
            chrom=AUTOSOME_LIST
        )
    output:
        metrics="results/sample_preqc/sample_metrics_by_chromosome.tsv",
        summary="results/sample_preqc/sample_chromosome_consistency.tsv"
    params:
        stats_dir="results/sample_preqc/by_chromosome",
        robust_z_threshold=config["sample_preqc"]["robust_z_threshold"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail

        python workflow/exploratory/sample_preqc_by_chromosome.py \
            --stats-dir {params.stats_dir} \
            --metrics-out {output.metrics} \
            --summary-out {output.summary} \
            --robust-z-threshold {params.robust_z_threshold}
        """
