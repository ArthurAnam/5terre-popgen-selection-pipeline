# ============================================================
# Cinque Terre Population Genomics and Selection Pipeline
# Main Snakefile skeleton
# ============================================================

configfile: "config/config.yaml"


# ============================================================
# Rule: all
# Final targets for the non-optional core workflow
# ============================================================

rule all:
    input:
        "results/harmonization/harmonization_done.txt",
        "results/qc/qc_done.txt",
        "results/relatedness/relatedness_done.txt",
        "results/pca/pca_done.txt",
        "results/roh/roh_done.txt",
        "results/ld_decay/ld_done.txt",
        "results/phasing/phasing_done.txt",
        "results/lassi/lassi_done.txt"


# ============================================================
# Dataset harmonization
# ============================================================

rule harmonization:
    output:
        "results/harmonization/harmonization_done.txt"
    shell:
        """
        echo "Dataset harmonization completed" > {output}
        """


# ============================================================
# Population-genomics QC
# ============================================================

rule qc:
    input:
        "results/harmonization/harmonization_done.txt"
    output:
        "results/qc/qc_done.txt"
    shell:
        """
        echo "Population-genomics QC completed" > {output}
        """


# ============================================================
# IBS / relatedness exploration
# ============================================================

rule relatedness:
    input:
        "results/qc/qc_done.txt"
    output:
        "results/relatedness/relatedness_done.txt"
    shell:
        """
        echo "IBS and relatedness exploration completed" > {output}
        """


# ============================================================
# PCA
# ============================================================

rule pca:
    input:
        qc="results/qc/qc_done.txt",
        relatedness="results/relatedness/relatedness_done.txt"
    output:
        "results/pca/pca_done.txt"
    shell:
        """
        echo "PCA completed" > {output}
        """


# ============================================================
# Runs of Homozygosity
# ============================================================

rule roh:
    input:
        qc="results/qc/qc_done.txt",
        relatedness="results/relatedness/relatedness_done.txt"
    output:
        "results/roh/roh_done.txt"
    shell:
        """
        echo "ROH analyses completed" > {output}
        """


# ============================================================
# LD decay
# ============================================================

rule ld_decay:
    input:
        "results/qc/qc_done.txt"
    output:
        "results/ld_decay/ld_done.txt"
    shell:
        """
        echo "LD decay estimation completed" > {output}
        """


# ============================================================
# Phasing
# ============================================================

rule phasing:
    input:
        "results/qc/qc_done.txt"
    output:
        "results/phasing/phasing_done.txt"
    shell:
        """
        echo "Phasing completed" > {output}
        """


# ============================================================
# LASSI selection scan
# ============================================================

rule lassi:
    input:
        ld_decay="results/ld_decay/ld_done.txt",
        phasing="results/phasing/phasing_done.txt"
    output:
        "results/lassi/lassi_done.txt"
    shell:
        """
        echo "LASSI scan completed" > {output}
        """
