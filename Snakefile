# ============================================================
# 5Terre Population Genomics Pipeline
# Main Snakefile (skeleton)
# ============================================================

configfile: "config/config.yaml"


# ============================================================
# Rule: all (final targets)
# ============================================================

rule all:
    input:
        "results/qc/qc_done.txt",
        "results/relatedness/master_unrelated.txt",
        "results/pca/pca_done.txt",
        "results/roh/roh_done.txt",
        "results/ld_decay/ld_done.txt",
        "results/phasing/phasing_done.txt",
        "results/lassi/lassi_done.txt"


# ============================================================
# QC
# ============================================================

rule qc:
    output:
        "results/qc/qc_done.txt"
    shell:
        """
        echo "QC completed" > {output}
        """


# ============================================================
# Relatedness
# ============================================================

rule relatedness:
    input:
        "results/qc/qc_done.txt"
    output:
        "results/relatedness/master_unrelated.txt"
    shell:
        """
        echo "Relatedness filtering completed" > {output}
        """


# ============================================================
# PCA
# ============================================================

rule pca:
    input:
        "results/relatedness/master_unrelated.txt"
    output:
        "results/pca/pca_done.txt"
    shell:
        """
        echo "PCA completed" > {output}
        """


# ============================================================
# ROH
# ============================================================

rule roh:
    input:
        "results/relatedness/master_unrelated.txt"
    output:
        "results/roh/roh_done.txt"
    shell:
        """
        echo "ROH completed" > {output}
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
        echo "LD decay completed" > {output}
        """


# ============================================================
# Phasing
# ============================================================

rule phasing:
    input:
        "results/ld_decay/ld_done.txt"
    output:
        "results/phasing/phasing_done.txt"
    shell:
        """
        echo "Phasing completed" > {output}
        """


# ============================================================
# LASSI
# ============================================================

rule lassi:
    input:
        "results/phasing/phasing_done.txt"
    output:
        "results/lassi/lassi_done.txt"
    shell:
        """
        echo "LASSI completed" > {output}
        """
