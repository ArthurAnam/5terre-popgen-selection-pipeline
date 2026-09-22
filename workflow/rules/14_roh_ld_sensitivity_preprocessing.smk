# ============================================================
# ROH LD-sensitivity preprocessing
# ============================================================
# Build the common MAF>=0.05 PLINK1 binary panel, apply the exact light
# VIF-pruning procedure used by Howrigan et al. (2011), and audit the physical
# density of the retained marker set before any final ROH parameters are frozen.

ROH_LD_SENS_DIR = "results/roh/ld_sensitivity"
ROH_COMMON_BED_PREFIX = ROH_LD_SENS_DIR + "/ct_1kg_eur.maf0.05"
ROH_LIGHT_VIF_PREFIX = ROH_LD_SENS_DIR + "/ct_1kg_eur.maf0.05.light_vif_50_5_10"
ROH_LIGHT_VIF_ORDERED_SNPLIST = ROH_LD_SENS_DIR + "/light_vif_50_5_10.ordered.snplist"


rule roh_common_maf_bed:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        snplist=ROH_MAF_SNPLIST
    output:
        bed=ROH_COMMON_BED_PREFIX + ".bed",
        bim=ROH_COMMON_BED_PREFIX + ".bim",
        fam=ROH_COMMON_BED_PREFIX + ".fam"
    log:
        "logs/roh/ld_sensitivity/common_maf_bed.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_LD_SENS_DIR} logs/roh/ld_sensitivity

        tmp_prefix="{ROH_LD_SENS_DIR}/.common_maf_bed.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.snplist} \
            --make-bed \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.bed" {output.bed}
        mv "$tmp_prefix.bim" {output.bim}
        mv "$tmp_prefix.fam" {output.fam}
        """


rule roh_light_vif_prune:
    input:
        bed=ROH_COMMON_BED_PREFIX + ".bed",
        bim=ROH_COMMON_BED_PREFIX + ".bim",
        fam=ROH_COMMON_BED_PREFIX + ".fam"
    output:
        prune_in=ROH_LIGHT_VIF_PREFIX + ".prune.in",
        prune_out=ROH_LIGHT_VIF_PREFIX + ".prune.out"
    params:
        window=lambda wildcards: config["population_structure"]["roh"]["ld_sensitivity"]["light_vif"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["roh"]["ld_sensitivity"]["light_vif"]["step_snps"],
        vif=lambda wildcards: config["population_structure"]["roh"]["ld_sensitivity"]["light_vif"]["vif_threshold"]
    log:
        "logs/roh/ld_sensitivity/light_vif_prune.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        tmp_prefix="{ROH_LD_SENS_DIR}/.light_vif.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink \
            --bfile {ROH_COMMON_BED_PREFIX} \
            --indep {params.window} {params.step} {params.vif} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule roh_light_vif_ordered_snplist:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=ROH_LIGHT_VIF_PREFIX + ".prune.in"
    output:
        ROH_LIGHT_VIF_ORDERED_SNPLIST
    log:
        "logs/roh/ld_sensitivity/light_vif_ordered_snplist.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        tmp_prefix="{ROH_LD_SENS_DIR}/.light_vif_ordered.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.prune_in} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output}
        """


rule audit_roh_light_vif_density:
    input:
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        snplist=ROH_LIGHT_VIF_ORDERED_SNPLIST,
        script="workflow/scripts/summarize_roh_marker_density.py"
    output:
        by_chr=ROH_LD_SENS_DIR + "/light_vif_marker_density_by_chromosome.tsv",
        gaps=ROH_LD_SENS_DIR + "/light_vif_intermarker_gap_thresholds.tsv",
        summary=ROH_LD_SENS_DIR + "/light_vif_marker_density_summary.tsv",
        readme=ROH_LD_SENS_DIR + "/README.txt"
    params:
        maf=lambda wildcards: config["population_structure"]["roh"]["maf_threshold"],
        expected_samples=lambda wildcards: config["population_structure"]["roh"]["expected_joint_samples"],
        gap_thresholds=lambda wildcards: ",".join(
            str(x) for x in config["population_structure"]["roh"]["marker_density_audit"]["gap_thresholds_kb"]
        )
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        n_variants=$(wc -l < {input.snplist})

        python {input.script} \
            --pvar {input.pvar} \
            --snplist {input.snplist} \
            --maf-threshold {params.maf} \
            --expected-variants "$n_variants" \
            --expected-samples {params.expected_samples} \
            --panel-label joint_harmonized_ct_plus_1kg_eur_maf0.05_light_vif \
            --ld-pruning-label PLINK1.9_indep_50_5_10_VIF \
            --gap-thresholds-kb {params.gap_thresholds} \
            --by-chromosome-out {output.by_chr} \
            --gap-thresholds-out {output.gaps} \
            --summary-out {output.summary} \
            --readme-out {output.readme}
        """


rule audit_roh_ld_sensitivity_preprocessing:
    input:
        ROH_LIGHT_VIF_PREFIX + ".prune.in",
        ROH_LIGHT_VIF_PREFIX + ".prune.out",
        ROH_LIGHT_VIF_ORDERED_SNPLIST,
        ROH_LD_SENS_DIR + "/light_vif_marker_density_by_chromosome.tsv",
        ROH_LD_SENS_DIR + "/light_vif_intermarker_gap_thresholds.tsv",
        ROH_LD_SENS_DIR + "/light_vif_marker_density_summary.tsv",
        ROH_LD_SENS_DIR + "/README.txt"
