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
        final_psam=PCA_PRUNED_PREFIX + ".psam"
    output:
        summary="results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
        maf_scope=lambda wildcards: config["population_structure"]["pca"]["maf_scope"],
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"],
        expected_samples=lambda wildcards: config["population_structure"]["pca"]["expected_joint_samples"],
        expected_variants=lambda wildcards: config["population_structure"]["pca"]["expected_harmonized_variants"]
    shell:
        r"""
        set -euo pipefail

        input_variants=$(grep -vc '^#' {input.input_pvar})
        input_samples=$(grep -vc '^#' {input.input_psam})
        maf_variants=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS
    output:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"],
        expected_samples=lambda wildcards: config["population_structure"]["pca"]["expected_joint_samples"],
        expected_variants=lambda wildcards: config["population_structure"]["pca"]["expected_harmonized_variants"]
    shell:
        r"""
        set -euo pipefail

        input_variants=$(grep -vc '^#' {input.input_pvar})
        input_samples=$(grep -vc '^#' {input.input_psam})
        maf_variants=$(grep -cve '^[[:space:]]*rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.maf_snplist})
        prune_in=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_in})
        prune_out=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_out})
        final_variants=$(grep -vc '^#' {input.final_pvar})
        final_samples=$(grep -vc '^#' {input.final_psam})

        [ "$input_variants" -eq "{params.expected_variants}" ] || {{
            echo "ERROR: expected {params.expected_variants} harmonized variants, found $input_variants" >&2
            exit 1
        }}
        [ "$input_samples" -eq "{params.expected_samples}" ] || {{
            echo "ERROR: expected {params.expected_samples} harmonized samples, found $input_samples" >&2
            exit 1
        }}
        [ $((prune_in + prune_out)) -eq "$maf_variants" ] || {{
            echo "ERROR: prune.in + prune.out does not equal MAF-passing marker count" >&2
            exit 1
        }}
        [ "$final_variants" -eq "$prune_in" ] || {{
            echo "ERROR: final pruned PVAR count does not equal prune.in count" >&2
            exit 1
        }}
        [ "$final_samples" -eq "$input_samples" ] || {{
            echo "ERROR: PCA preprocessing changed sample count" >&2
            exit 1
        }}

        removed_maf=$((input_variants - maf_variants))

        {{
            printf 'metric\tvalue\n'
            printf 'harmonized_input_variants\t%s\n' "$input_variants"
            printf 'harmonized_input_samples\t%s\n' "$input_samples"
            printf 'maf_threshold\t%s\n' "{params.maf}"
            printf 'maf_scope\t%s\n' "{params.maf_scope}"
            printf 'variants_after_maf\t%s\n' "$maf_variants"
            printf 'variants_removed_by_maf\t%s\n' "$removed_maf"
            printf 'ld_window_snps\t%s\n' "{params.window}"
            printf 'ld_step_snps\t%s\n' "{params.step}"
            printf 'ld_r2_threshold\t%s\n' "{params.r2}"
            printf 'plink2_indep_order\t%s\n' "{params.indep_order}"
            printf 'variants_removed_by_ld_pruning\t%s\n' "$prune_out"
            printf 'variants_after_ld_pruning\t%s\n' "$prune_in"
            printf 'pca_preprocessed_samples\t%s\n' "$final_samples"
            printf 'pca_preprocessed_variants\t%s\n' "$final_variants"
        }} > {output.summary}
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.maf_snplist})
        masked_variants=$(grep -cve '^[[:space:]]*rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.maf_snplist})
        prune_in=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_in})
        prune_out=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_out})
        final_variants=$(grep -vc '^#' {input.final_pvar})
        final_samples=$(grep -vc '^#' {input.final_psam})

        [ "$input_variants" -eq "{params.expected_variants}" ] || {{
            echo "ERROR: expected {params.expected_variants} harmonized variants, found $input_variants" >&2
            exit 1
        }}
        [ "$input_samples" -eq "{params.expected_samples}" ] || {{
            echo "ERROR: expected {params.expected_samples} harmonized samples, found $input_samples" >&2
            exit 1
        }}
        [ $((prune_in + prune_out)) -eq "$maf_variants" ] || {{
            echo "ERROR: prune.in + prune.out does not equal MAF-passing marker count" >&2
            exit 1
        }}
        [ "$final_variants" -eq "$prune_in" ] || {{
            echo "ERROR: final pruned PVAR count does not equal prune.in count" >&2
            exit 1
        }}
        [ "$final_samples" -eq "$input_samples" ] || {{
            echo "ERROR: PCA preprocessing changed sample count" >&2
            exit 1
        }}

        removed_maf=$((input_variants - maf_variants))

        {{
            printf 'metric\tvalue\n'
            printf 'harmonized_input_variants\t%s\n' "$input_variants"
            printf 'harmonized_input_samples\t%s\n' "$input_samples"
            printf 'maf_threshold\t%s\n' "{params.maf}"
            printf 'maf_scope\t%s\n' "{params.maf_scope}"
            printf 'variants_after_maf\t%s\n' "$maf_variants"
            printf 'variants_removed_by_maf\t%s\n' "$removed_maf"
            printf 'ld_window_snps\t%s\n' "{params.window}"
            printf 'ld_step_snps\t%s\n' "{params.step}"
            printf 'ld_r2_threshold\t%s\n' "{params.r2}"
            printf 'plink2_indep_order\t%s\n' "{params.indep_order}"
            printf 'variants_removed_by_ld_pruning\t%s\n' "$prune_out"
            printf 'variants_after_ld_pruning\t%s\n' "$prune_in"
            printf 'pca_preprocessed_samples\t%s\n' "$final_samples"
            printf 'pca_preprocessed_variants\t%s\n' "$final_variants"
        }} > {output.summary}
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.masked_snplist})
        unmasked_prune_in=$(grep -cve '^[[:space:]]*rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.maf_snplist})
        prune_in=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_in})
        prune_out=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_out})
        final_variants=$(grep -vc '^#' {input.final_pvar})
        final_samples=$(grep -vc '^#' {input.final_psam})

        [ "$input_variants" -eq "{params.expected_variants}" ] || {{
            echo "ERROR: expected {params.expected_variants} harmonized variants, found $input_variants" >&2
            exit 1
        }}
        [ "$input_samples" -eq "{params.expected_samples}" ] || {{
            echo "ERROR: expected {params.expected_samples} harmonized samples, found $input_samples" >&2
            exit 1
        }}
        [ $((prune_in + prune_out)) -eq "$maf_variants" ] || {{
            echo "ERROR: prune.in + prune.out does not equal MAF-passing marker count" >&2
            exit 1
        }}
        [ "$final_variants" -eq "$prune_in" ] || {{
            echo "ERROR: final pruned PVAR count does not equal prune.in count" >&2
            exit 1
        }}
        [ "$final_samples" -eq "$input_samples" ] || {{
            echo "ERROR: PCA preprocessing changed sample count" >&2
            exit 1
        }}

        removed_maf=$((input_variants - maf_variants))

        {{
            printf 'metric\tvalue\n'
            printf 'harmonized_input_variants\t%s\n' "$input_variants"
            printf 'harmonized_input_samples\t%s\n' "$input_samples"
            printf 'maf_threshold\t%s\n' "{params.maf}"
            printf 'maf_scope\t%s\n' "{params.maf_scope}"
            printf 'variants_after_maf\t%s\n' "$maf_variants"
            printf 'variants_removed_by_maf\t%s\n' "$removed_maf"
            printf 'ld_window_snps\t%s\n' "{params.window}"
            printf 'ld_step_snps\t%s\n' "{params.step}"
            printf 'ld_r2_threshold\t%s\n' "{params.r2}"
            printf 'plink2_indep_order\t%s\n' "{params.indep_order}"
            printf 'variants_removed_by_ld_pruning\t%s\n' "$prune_out"
            printf 'variants_after_ld_pruning\t%s\n' "$prune_in"
            printf 'pca_preprocessed_samples\t%s\n' "$final_samples"
            printf 'pca_preprocessed_variants\t%s\n' "$final_variants"
        }} > {output.summary}
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.unmasked_prune_in})
        masked_prune_in=$(grep -cve '^[[:space:]]*rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.maf_snplist})
        prune_in=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_in})
        prune_out=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_out})
        final_variants=$(grep -vc '^#' {input.final_pvar})
        final_samples=$(grep -vc '^#' {input.final_psam})

        [ "$input_variants" -eq "{params.expected_variants}" ] || {{
            echo "ERROR: expected {params.expected_variants} harmonized variants, found $input_variants" >&2
            exit 1
        }}
        [ "$input_samples" -eq "{params.expected_samples}" ] || {{
            echo "ERROR: expected {params.expected_samples} harmonized samples, found $input_samples" >&2
            exit 1
        }}
        [ $((prune_in + prune_out)) -eq "$maf_variants" ] || {{
            echo "ERROR: prune.in + prune.out does not equal MAF-passing marker count" >&2
            exit 1
        }}
        [ "$final_variants" -eq "$prune_in" ] || {{
            echo "ERROR: final pruned PVAR count does not equal prune.in count" >&2
            exit 1
        }}
        [ "$final_samples" -eq "$input_samples" ] || {{
            echo "ERROR: PCA preprocessing changed sample count" >&2
            exit 1
        }}

        removed_maf=$((input_variants - maf_variants))

        {{
            printf 'metric\tvalue\n'
            printf 'harmonized_input_variants\t%s\n' "$input_variants"
            printf 'harmonized_input_samples\t%s\n' "$input_samples"
            printf 'maf_threshold\t%s\n' "{params.maf}"
            printf 'maf_scope\t%s\n' "{params.maf_scope}"
            printf 'variants_after_maf\t%s\n' "$maf_variants"
            printf 'variants_removed_by_maf\t%s\n' "$removed_maf"
            printf 'ld_window_snps\t%s\n' "{params.window}"
            printf 'ld_step_snps\t%s\n' "{params.step}"
            printf 'ld_r2_threshold\t%s\n' "{params.r2}"
            printf 'plink2_indep_order\t%s\n' "{params.indep_order}"
            printf 'variants_removed_by_ld_pruning\t%s\n' "$prune_out"
            printf 'variants_after_ld_pruning\t%s\n' "$prune_in"
            printf 'pca_preprocessed_samples\t%s\n' "$final_samples"
            printf 'pca_preprocessed_variants\t%s\n' "$final_variants"
        }} > {output.summary}
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.masked_prune_in})
        masked_prune_out=$(grep -cve '^[[:space:]]*rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.maf_snplist})
        prune_in=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_in})
        prune_out=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_out})
        final_variants=$(grep -vc '^#' {input.final_pvar})
        final_samples=$(grep -vc '^#' {input.final_psam})

        [ "$input_variants" -eq "{params.expected_variants}" ] || {{
            echo "ERROR: expected {params.expected_variants} harmonized variants, found $input_variants" >&2
            exit 1
        }}
        [ "$input_samples" -eq "{params.expected_samples}" ] || {{
            echo "ERROR: expected {params.expected_samples} harmonized samples, found $input_samples" >&2
            exit 1
        }}
        [ $((prune_in + prune_out)) -eq "$maf_variants" ] || {{
            echo "ERROR: prune.in + prune.out does not equal MAF-passing marker count" >&2
            exit 1
        }}
        [ "$final_variants" -eq "$prune_in" ] || {{
            echo "ERROR: final pruned PVAR count does not equal prune.in count" >&2
            exit 1
        }}
        [ "$final_samples" -eq "$input_samples" ] || {{
            echo "ERROR: PCA preprocessing changed sample count" >&2
            exit 1
        }}

        removed_maf=$((input_variants - maf_variants))

        {{
            printf 'metric\tvalue\n'
            printf 'harmonized_input_variants\t%s\n' "$input_variants"
            printf 'harmonized_input_samples\t%s\n' "$input_samples"
            printf 'maf_threshold\t%s\n' "{params.maf}"
            printf 'maf_scope\t%s\n' "{params.maf_scope}"
            printf 'variants_after_maf\t%s\n' "$maf_variants"
            printf 'variants_removed_by_maf\t%s\n' "$removed_maf"
            printf 'ld_window_snps\t%s\n' "{params.window}"
            printf 'ld_step_snps\t%s\n' "{params.step}"
            printf 'ld_r2_threshold\t%s\n' "{params.r2}"
            printf 'plink2_indep_order\t%s\n' "{params.indep_order}"
            printf 'variants_removed_by_ld_pruning\t%s\n' "$prune_out"
            printf 'variants_after_ld_pruning\t%s\n' "$prune_in"
            printf 'pca_preprocessed_samples\t%s\n' "$final_samples"
            printf 'pca_preprocessed_variants\t%s\n' "$final_variants"
        }} > {output.summary}
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.masked_prune_out})
        final_variants=$(grep -vc '^#' {input.final_pvar})
        final_samples=$(grep -vc '^#' {input.final_psam})
        region_count=$(grep -cvE '^[[:space:]]*($|#)' {input.regions})

        [ "$input_variants" -eq "{params.expected_variants}" ] || {{
            echo "ERROR: expected {params.expected_variants} harmonized variants, found $input_variants" >&2
            exit 1
        }}
        [ "$input_samples" -eq "{params.expected_samples}" ] || {{
            echo "ERROR: expected {params.expected_samples} harmonized samples, found $input_samples" >&2
            exit 1
        }}
        [ "$masked_variants" -le "$maf_variants" ] || {{
            echo "ERROR: high-LD masking increased marker count" >&2
            exit 1
        }}
        [ $((masked_prune_in + masked_prune_out)) -eq "$masked_variants" ] || {{
            echo "ERROR: masked prune.in + prune.out does not equal masked marker count" >&2
            exit 1
        }}
        [ "$final_variants" -eq "$masked_prune_in" ] || {{
            echo "ERROR: masked final PVAR count does not equal masked prune.in count" >&2
            exit 1
        }}
        [ "$final_samples" -eq "$input_samples" ] || {{
            echo "ERROR: high-LD sensitivity changed sample count" >&2
            exit 1
        }}

        removed_highld=$((maf_variants - masked_variants))
        masked_minus_unmasked=$((masked_prune_in - unmasked_prune_in))

        {{
            printf 'metric\tvalue\n'
            printf 'harmonized_input_variants\t%s\n' "$input_variants"
            printf 'harmonized_input_samples\t%s\n' "$input_samples"
            printf 'maf_threshold\t%s\n' "{params.maf}"
            printf 'variants_after_maf\t%s\n' "$maf_variants"
            printf 'high_ld_region_count\t%s\n' "$region_count"
            printf 'variants_removed_by_high_ld_mask\t%s\n' "$removed_highld"
            printf 'variants_entering_masked_ld_pruning\t%s\n' "$masked_variants"
            printf 'ld_window_snps\t%s\n' "{params.window}"
            printf 'ld_step_snps\t%s\n' "{params.step}"
            printf 'ld_r2_threshold\t%s\n' "{params.r2}"
            printf 'plink2_indep_order\t%s\n' "{params.indep_order}"
            printf 'variants_removed_by_masked_ld_pruning\t%s\n' "$masked_prune_out"
            printf 'masked_variants_after_ld_pruning\t%s\n' "$masked_prune_in"
            printf 'unmasked_variants_after_ld_pruning\t%s\n' "$unmasked_prune_in"
            printf 'masked_minus_unmasked_final_variants\t%s\n' "$masked_minus_unmasked"
            printf 'pca_preprocessed_samples_masked\t%s\n' "$final_samples"
            printf 'pca_preprocessed_variants_masked\t%s\n' "$final_variants"
        }} > {output.summary}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.maf_snplist})
        prune_in=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_in})
        prune_out=$(grep -cve '^[[:space:]]*rule write_joint_pca_preprocessing_readme:
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
 {input.prune_out})
        final_variants=$(grep -vc '^#' {input.final_pvar})
        final_samples=$(grep -vc '^#' {input.final_psam})

        [ "$input_variants" -eq "{params.expected_variants}" ] || {{
            echo "ERROR: expected {params.expected_variants} harmonized variants, found $input_variants" >&2
            exit 1
        }}
        [ "$input_samples" -eq "{params.expected_samples}" ] || {{
            echo "ERROR: expected {params.expected_samples} harmonized samples, found $input_samples" >&2
            exit 1
        }}
        [ $((prune_in + prune_out)) -eq "$maf_variants" ] || {{
            echo "ERROR: prune.in + prune.out does not equal MAF-passing marker count" >&2
            exit 1
        }}
        [ "$final_variants" -eq "$prune_in" ] || {{
            echo "ERROR: final pruned PVAR count does not equal prune.in count" >&2
            exit 1
        }}
        [ "$final_samples" -eq "$input_samples" ] || {{
            echo "ERROR: PCA preprocessing changed sample count" >&2
            exit 1
        }}

        removed_maf=$((input_variants - maf_variants))

        {{
            printf 'metric\tvalue\n'
            printf 'harmonized_input_variants\t%s\n' "$input_variants"
            printf 'harmonized_input_samples\t%s\n' "$input_samples"
            printf 'maf_threshold\t%s\n' "{params.maf}"
            printf 'maf_scope\t%s\n' "{params.maf_scope}"
            printf 'variants_after_maf\t%s\n' "$maf_variants"
            printf 'variants_removed_by_maf\t%s\n' "$removed_maf"
            printf 'ld_window_snps\t%s\n' "{params.window}"
            printf 'ld_step_snps\t%s\n' "{params.step}"
            printf 'ld_r2_threshold\t%s\n' "{params.r2}"
            printf 'plink2_indep_order\t%s\n' "{params.indep_order}"
            printf 'variants_removed_by_ld_pruning\t%s\n' "$prune_out"
            printf 'variants_after_ld_pruning\t%s\n' "$prune_in"
            printf 'pca_preprocessed_samples\t%s\n' "$final_samples"
            printf 'pca_preprocessed_variants\t%s\n' "$final_variants"
        }} > {output.summary}
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


# ============================================================
# Sensitivity branch: predefined high/long-range-LD mask
# ============================================================
#
# This branch deliberately keeps the already generated unmasked panel.
# It applies a predefined GRCh37 high-LD mask after the joint MAF filter,
# then reruns the same LD-pruning parameters.  The two PCA marker panels can
# therefore be compared without changing any other preprocessing choice.

PCA_HIGHLD_REGIONS = "config/pca_high_ld_regions_grch37.bed1"
PCA_HIGHLD_SNPLIST = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.snplist"
PCA_HIGHLD_PRUNE_IN = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.in"
PCA_HIGHLD_PRUNE_OUT = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.prune.out"
PCA_HIGHLD_PRUNED_PREFIX = PCA_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05.highld_masked.ldpruned"


rule pca_joint_highld_mask_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        regions=PCA_HIGHLD_REGIONS
    output:
        snplist=PCA_HIGHLD_SNPLIST
    log:
        "logs/population_structure/pca/highld_mask.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_mask.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.maf_snplist} \
            --exclude bed1 {input.regions} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output.snplist}
        """


rule pca_joint_highld_ld_pruning:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        masked_snplist=PCA_HIGHLD_SNPLIST
    output:
        prune_in=PCA_HIGHLD_PRUNE_IN,
        prune_out=PCA_HIGHLD_PRUNE_OUT
    params:
        window=lambda wildcards: config["population_structure"]["pca"]["pruning"]["window_snps"],
        step=lambda wildcards: config["population_structure"]["pca"]["pruning"]["step_snps"],
        r2=lambda wildcards: config["population_structure"]["pca"]["pruning"]["r2_threshold"],
        indep_order=lambda wildcards: config["population_structure"]["pca"]["pruning"]["indep_order"]
    log:
        "logs/population_structure/pca/highld_masked_ld_pruning.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.highld_masked_ld_pruning.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.masked_snplist} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --indep-order {params.indep_order} \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.prune.in" {output.prune_in}
        mv "$tmp_prefix.prune.out" {output.prune_out}
        """


rule make_joint_pca_highld_masked_pgen:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        prune_in=PCA_HIGHLD_PRUNE_IN
    output:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam"
    log:
        "logs/population_structure/pca/make_highld_masked_pgen.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_PREPROCESSING_DIR} logs/population_structure/pca

        tmp_prefix="{PCA_PREPROCESSING_DIR}/.make_highld_masked_pgen.$$"
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


rule summarize_joint_pca_highld_sensitivity:
    input:
        input_pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        input_psam=PCA_HARMONIZED_PREFIX + ".psam",
        maf_snplist=PCA_MAF_SNPLIST,
        masked_snplist=PCA_HIGHLD_SNPLIST,
        unmasked_prune_in=PCA_PRUNE_IN,
        masked_prune_in=PCA_HIGHLD_PRUNE_IN,
        masked_prune_out=PCA_HIGHLD_PRUNE_OUT,
        final_pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        final_psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        regions=PCA_HIGHLD_REGIONS,
        script="workflow/scripts/summarize_pca_highld_sensitivity.py"
    output:
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv"
    params:
        maf=lambda wildcards: config["population_structure"]["pca"]["maf_threshold"],
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
            --masked-snplist {input.masked_snplist} \
            --unmasked-prune-in {input.unmasked_prune_in} \
            --masked-prune-in {input.masked_prune_in} \
            --masked-prune-out {input.masked_prune_out} \
            --final-pvar {input.final_pvar} \
            --final-psam {input.final_psam} \
            --regions {input.regions} \
            --maf-threshold {params.maf} \
            --window-snps {params.window} \
            --step-snps {params.step} \
            --r2-threshold {params.r2} \
            --indep-order {params.indep_order} \
            --expected-samples {params.expected_samples} \
            --expected-input-variants {params.expected_variants} \
            --summary-out {output}
        """


rule write_joint_pca_highld_sensitivity_readme:
    input:
        summary="results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        regions=PCA_HIGHLD_REGIONS
    output:
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Joint PCA high-LD sensitivity branch
====================================

Purpose
-------
The unmasked PCA-preprocessing panel is retained unchanged.  This sensitivity
branch starts from the same 549-person harmonized dataset and the same joint
MAF >= 0.05 SNP list, then excludes a predefined GRCh37 high/long-range-LD
region catalog before rerunning exactly the same LD-pruning parameters.

The masked and unmasked analyses therefore differ only by the predefined
region mask.  They are not alternative data harmonizations and no result-driven
regions are added.

Region catalog
--------------
The interval file is:
    config/pca_high_ld_regions_grch37.bed1

It is the GRCh37 interval list distributed with Grinde et al. 2024 after their
literature review of high, long-range, or otherwise unusual LD regions.  This
choice avoids directly copying the coarse Price et al. 2008 Mb intervals across
genome builds.  See config/pca_high_ld_regions_grch37.README.md and the project
methodological references for provenance.

Filtering order
---------------
  1. Joint MAF >= 0.05 on all 549 individuals (shared with unmasked branch).
  2. Exclude predefined GRCh37 high-LD intervals.
  3. LD pruning with --indep-pairwise 50 5 0.2 --indep-order 1.

The final choice of which version is emphasized in the paper is intentionally
left pending until the masked and unmasked PCA results have been compared.
EOF
        """


rule prepare_joint_pca_sensitivity:
    input:
        "results/population_structure/pca/preprocessing/pca_preprocessing_summary.tsv",
        PCA_PRUNED_PREFIX + ".pgen",
        PCA_PRUNED_PREFIX + ".pvar",
        PCA_PRUNED_PREFIX + ".psam",
        PCA_HIGHLD_SNPLIST,
        PCA_HIGHLD_PRUNE_IN,
        PCA_HIGHLD_PRUNE_OUT,
        PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        "results/population_structure/pca/preprocessing/pca_highld_sensitivity_summary.tsv",
        "results/population_structure/pca/preprocessing/HIGH_LD_SENSITIVITY_README.txt"
