# ============================================================
# Joint PCA preprocessing: MAF filter + LD pruning
# ============================================================
#
# Input is the already harmonized Cinque Terre + 1000G EUR PGEN dataset.
# The PCA marker panel is defined once on the joint 549-sample dataset:
#   1) MAF >= configured threshold
#   2) LD pruning on the MAF-filtered marker set
#
# No separate target/reference pruning lists are generated.

PCA_HARMONIZED_PREFIX = "results/population_structure/harmonization/ct_1kg_eur.harmonized"
PCA_PREPROCESSING_DIR = "results/population_structure/pca/preprocessing"
PCA_MAF_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.snplist"
PCA_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.prune.in"
PCA_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.prune.out"
PCA_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.ldpruned"


rule pca_joint_maf_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam"
    output:
        snplist=PCA_MAF_SNPLIST
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"]
    log:
        "logs/population_structure/pca/maf_filter.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.maf_filter.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --maf {params.maf} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST
    output:
        prune_in=PCA_PRUNE_IN,
        prune_out=PCA_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_pruned_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_PRUNE_IN
    output:
        pgen=PCA_PRUNED_PREFIX + ".pgen",
        pvar=PCA_PRUNED_PREFIX + ".pvar",
        psam=PCA_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_pruned_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_pruned_pgen.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.prune_in} \
            --make-pgen \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.pgen" {output.pgen}
        mv "$tmp_prefix.pvar" {output.pvar}
        mv "$tmp_prefix.psam" {output.psam}
        """


rule summarize_joint_pca_preprocessing:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        prune_in=PCA_PRUNE_IN,
        prune_out=PCA_PRUNE_OUT,
        final_pvar=PCA_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_PRUNED_PREFIX + ".psam",
        script="workflow/scripts/summarize_pca_preprocessing.py"
    output:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
        maf_scope=lambda wildcards: config["population_structure"]["pca"]["maf_scope"],
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"],
        expected_samples=lambda wildcards: config["population_structure"]["pca"]["expected_joint_samples"],
        expected_variants=lambda wildcards: config["population_structure"]["pca"]["expected_harmonized_variants"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --input-pvar {input.input_pvar} \
            --input-psam {input.input_psam} \
            --maf-snplist {input.maf_snplist} \
            --prune-in {input.prune_in} \
            --prune-out {input.prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --maf-threshold {params.maf} \
            --maf-scope {params.maf_scope} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_preprocessing_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv"
    output:
        "results/population_structure/pca/preprocessing/README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA preprocessing: Cinque Terre + 1000 Genomes EUR
=========================================================

Input
-----
The input is the persistent harmonized dataset containing 46 Cinque Terre and
503 1000 Genomes EUR individuals (549 total), restricted to exact
CHR:POS:REF:ALT matches.

Marker filtering
----------------
The PCA marker panel is defined once on the JOINT harmonized dataset. There are
not separate Cinque Terre and 1000G marker-pruning branches.

Filtering order:
  1. MAF >= 0.05 on all 549 individuals.
  2. LD pruning on the MAF-filtered marker set with PLINK2:
       --indep-pairwise 50 5 0.2
       --indep-order 1

The explicit --indep-order 1 setting is used to retain PLINK 1.x pruning-order
behavior instead of relying on the newer PLINK2 default.

Persistent outputs
------------------
The MAF-passing SNP list and the prune.in/prune.out lists are retained as audit
artifacts. The LD-pruned PGEN/PVAR/PSAM dataset is also persistent.

This rule set prepares the marker dataset only. smartpca itself is not run here;
the remaining smartpca-specific decisions (including number of PCs, explicit
long-range-LD masking, and outlier handling) are handled separately.
EOF
        """


rule prepare_joint_pca:
    input:
        PCA_MAF_SNPLIST,
        PCA_PRUNE_IN,
        PCA_PRUNE_OUT,
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        "results/population_structure/pca/preprocessing/README.txt"
