# ============================================================
# Individual pre-QC: PLINK2 heterozygosity confirmation
# ============================================================
# Heterozygosity is re-estimated on a common, LD-pruned autosomal marker set.
# This is diagnostic only. No sample is removed by these rules.

HET_CFG = config["sample_preqc"]["heterozygosity_confirmation"]
HET_LD = HET_CFG["ld_pruning"]


rule sample_preqc_plink2_marker_set:
    input:
        vcf=lambda wildcards: config["local_paths"]["target_vcf"],
        preqc="results/preqc/preqc_done.txt"
    output:
        pgen=temp("results/sample_preqc/plink2/common_autosomal.pgen"),
        pvar=temp("results/sample_preqc/plink2/common_autosomal.pvar"),
        psam=temp("results/sample_preqc/plink2/common_autosomal.psam")
    log:
        "logs/sample_preqc/plink2_marker_set.log"
    params:
        maf=HET_CFG["maf_threshold"],
        geno=HET_CFG["variant_missingness_threshold"],
        out="results/sample_preqc/plink2/common_autosomal"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/sample_preqc/plink2 logs/sample_preqc

        plink2 \
            --vcf {input.vcf} \
            --chr 1-22 \
            --snps-only just-acgt \
            --max-alleles 2 \
            --maf {params.maf} \
            --geno {params.geno} \
            --set-all-var-ids '@:#:$r:$a' \
            --make-pgen \
            --out {params.out} \
            > {log} 2>&1
        """


rule sample_preqc_plink2_ld_prune:
    input:
        pgen="results/sample_preqc/plink2/common_autosomal.pgen",
        pvar="results/sample_preqc/plink2/common_autosomal.pvar",
        psam="results/sample_preqc/plink2/common_autosomal.psam"
    output:
        prune_in="results/sample_preqc/plink2/common_autosomal.prune.in",
        prune_out="results/sample_preqc/plink2/common_autosomal.prune.out"
    log:
        "logs/sample_preqc/plink2_ld_prune.log"
    params:
        window=HET_LD["window_variants"],
        step=HET_LD["step_variants"],
        r2=HET_LD["r2_threshold"],
        pfile="results/sample_preqc/plink2/common_autosomal",
        out="results/sample_preqc/plink2/common_autosomal"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail

        plink2 \
            --pfile {params.pfile} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --out {params.out} \
            > {log} 2>&1
        """


rule sample_preqc_plink2_heterozygosity:
    input:
        pgen="results/sample_preqc/plink2/common_autosomal.pgen",
        pvar="results/sample_preqc/plink2/common_autosomal.pvar",
        psam="results/sample_preqc/plink2/common_autosomal.psam",
        prune_in="results/sample_preqc/plink2/common_autosomal.prune.in"
    output:
        het="results/sample_preqc/plink2/heterozygosity.het"
    log:
        "logs/sample_preqc/plink2_heterozygosity.log"
    params:
        pfile="results/sample_preqc/plink2/common_autosomal",
        out="results/sample_preqc/plink2/heterozygosity"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail

        plink2 \
            --pfile {params.pfile} \
            --extract {input.prune_in} \
            --het \
            --out {params.out} \
            > {log} 2>&1
        """
