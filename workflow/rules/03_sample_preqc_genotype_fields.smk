# ============================================================
# Individual pre-QC: genotype-field diagnostics
# ============================================================
# Uses a deterministic genome-wide subset of autosomal regions to compare
# GT/DP/GQ/AD distributions across samples. No external hard thresholds are
# applied and no sample is removed automatically.

GENOTYPE_QC_REGIONS = ",".join(
    f"{chrom}:{start}-{start + 999999}"
    for chrom in range(1, 23)
    for start in (20000000, 35000000, 45000000)
)


rule sample_preqc_genotype_fields:
    input:
        vcf=lambda wildcards: config["local_paths"]["target_vcf"],
        samples="results/preqc/samples.txt"
    output:
        metrics="results/sample_preqc/sample_genotype_field_diagnostics.tsv"
    log:
        "logs/sample_preqc/genotype_field_diagnostics.log"
    params:
        regions=GENOTYPE_QC_REGIONS
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/sample_preqc logs/sample_preqc

        bcftools query \
            -r {params.regions} \
            -f '%CHROM\t%POS[\t%GT\t%DP\t%GQ\t%AD]\n' \
            {input.vcf} 2> {log} \
        | python workflow/scripts/sample_preqc_genotype_fields.py \
            --samples {input.samples} \
            --metrics-out {output.metrics}
        """
