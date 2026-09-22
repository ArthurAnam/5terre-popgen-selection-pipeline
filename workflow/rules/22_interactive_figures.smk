# ============================================================
# Cross-analysis interactive figures and caption registry
# ============================================================

FIG_DIR = "results/figures"

rule build_interactive_analysis_figures:
    input:
        registry="config/figure_registry.tsv",
        king_pairs="results/relatedness/king/king_pairwise_kinship.tsv",
        pca_unmasked=PCA_SMARTPCA_DIR + "/unmasked/pca_coordinates.tsv",
        pca_unmasked_eigen=PCA_SMARTPCA_DIR + "/unmasked/pca_eigenvalues.tsv",
        pca_masked=PCA_SMARTPCA_DIR + "/highld_masked/pca_coordinates.tsv",
        pca_masked_eigen=PCA_SMARTPCA_DIR + "/highld_masked/pca_eigenvalues.tsv",
        roh_individual=ROH_PROD_DIR + "/individual_roh_froh.tsv",
        ld_maf005_bins=SEL_LD_DIR + "/summary/maf005.ld_decay_bins.tsv",
        ld_maf005_summary=SEL_LD_DIR + "/summary/maf005.ld_decay_summary.tsv",
        ld_maf001_bins=SEL_LD_DIR + "/summary/maf001.ld_decay_bins.tsv",
        ld_maf001_summary=SEL_LD_DIR + "/summary/maf001.ld_decay_summary.tsv",
        script="workflow/scripts/build_interactive_figures.py"
    output:
        king=FIG_DIR + "/king_kinship_vs_ibs0.html",
        pca12=FIG_DIR + "/pca_unmasked_pc1_pc2.html",
        pca23=FIG_DIR + "/pca_unmasked_pc2_pc3.html",
        pca12_masked=FIG_DIR + "/pca_highld_pc1_pc2.html",
        pca23_masked=FIG_DIR + "/pca_highld_pc2_pc3.html",
        froh=FIG_DIR + "/roh_froh_ge1_5mb.html",
        froh_long=FIG_DIR + "/roh_froh_ge5mb.html",
        roh_scatter=FIG_DIR + "/roh_nroh_vs_total.html",
        ld_primary=FIG_DIR + "/ld_decay_maf005.html",
        ld_sensitivity=FIG_DIR + "/ld_decay_maf001.html",
        index=FIG_DIR + "/index.html"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {FIG_DIR}
        python {input.script} \
            --registry {input.registry} \
            --king-pairs {input.king_pairs} \
            --pca-unmasked {input.pca_unmasked} \
            --pca-unmasked-eigen {input.pca_unmasked_eigen} \
            --pca-masked {input.pca_masked} \
            --pca-masked-eigen {input.pca_masked_eigen} \
            --roh-individual {input.roh_individual} \
            --ld-maf005-bins {input.ld_maf005_bins} \
            --ld-maf005-summary {input.ld_maf005_summary} \
            --ld-maf001-bins {input.ld_maf001_bins} \
            --ld-maf001-summary {input.ld_maf001_summary} \
            --out-dir {FIG_DIR} \
            --index-out {output.index}
        """

rule interactive_analysis_figures:
    input:
        rules.build_interactive_analysis_figures.output
