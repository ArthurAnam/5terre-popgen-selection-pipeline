# ============================================================
# Joint smartpca: unmasked and high-LD-masked sensitivity panels
# ============================================================

PCA_SMARTPCA_DIR = "results/population_structure/pca/smartpca"
PCA_SAMPLE_ANNOTATIONS = PCA_SMARTPCA_DIR + "/sample_population_annotations.tsv"
PCA_SAMPLE_COUNTS = PCA_SMARTPCA_DIR + "/sample_population_counts.tsv"

PCA_SMARTPCA_PANEL_PREFIXES = {
    "unmasked": PCA_PRUNED_PREFIX,
    "highld_masked": PCA_HIGHLD_PRUNED_PREFIX,
}


def pca_smartpca_panel_prefix(wildcards):
    try:
        return PCA_SMARTPCA_PANEL_PREFIXES[wildcards.pca_panel]
    except KeyError as exc:
        raise ValueError(f"Unknown PCA panel: {wildcards.pca_panel}") from exc


wildcard_constraints:
    pca_panel="unmasked|highld_masked"


rule prepare_pca_sample_annotations:
    input:
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        eur_metadata="results/reference/1kg_eur/metadata/eur_samples.tsv",
        script="workflow/scripts/prepare_pca_sample_annotations.py"
    output:
        annotations=PCA_SAMPLE_ANNOTATIONS,
        counts=PCA_SAMPLE_COUNTS
    params:
        expected_total=lambda wildcards: config["population_structure"]["pca"]["expected_joint_samples"],
        expected_eur=lambda wildcards: config["datasets"]["reference_1kg_eur"]["populations"]["expected_total_eur_samples"],
        expected_ct=46
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_SMARTPCA_DIR}
        python {input.script} \
            --psam {input.psam} \
            --eur-metadata {input.eur_metadata} \
            --expected-total {params.expected_total} \
            --expected-eur {params.expected_eur} \
            --expected-ct {params.expected_ct} \
            --annotations-out {output.annotations} \
            --counts-out {output.counts}
        """


rule make_smartpca_plink_input:
    input:
        pgen=lambda wildcards: pca_smartpca_panel_prefix(wildcards) + ".pgen",
        pvar=lambda wildcards: pca_smartpca_panel_prefix(wildcards) + ".pvar",
        psam=lambda wildcards: pca_smartpca_panel_prefix(wildcards) + ".psam",
        annotations=PCA_SAMPLE_ANNOTATIONS,
        script="workflow/scripts/label_smartpca_fam.py"
    output:
        bed=temp(PCA_SMARTPCA_DIR + "/{pca_panel}/input.bed"),
        bim=temp(PCA_SMARTPCA_DIR + "/{pca_panel}/input.bim"),
        fam=temp(PCA_SMARTPCA_DIR + "/{pca_panel}/input.fam")
    log:
        "logs/population_structure/pca/smartpca/{pca_panel}.plink2.log"
    params:
        prefix=pca_smartpca_panel_prefix
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_SMARTPCA_DIR}/{wildcards.pca_panel} logs/population_structure/pca/smartpca

        tmp_prefix="{PCA_SMARTPCA_DIR}/{wildcards.pca_panel}/.plink_input.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {params.prefix} \
            --make-bed \
            --out "$tmp_prefix" \
            > {log} 2>&1

        python {input.script} \
            --fam "$tmp_prefix.fam" \
            --annotations {input.annotations} \
            --output-fam {output.fam}

        mv "$tmp_prefix.bed" {output.bed}
        mv "$tmp_prefix.bim" {output.bim}
        """


rule run_joint_smartpca:
    input:
        bed=PCA_SMARTPCA_DIR + "/{pca_panel}/input.bed",
        bim=PCA_SMARTPCA_DIR + "/{pca_panel}/input.bim",
        fam=PCA_SMARTPCA_DIR + "/{pca_panel}/input.fam"
    output:
        par=PCA_SMARTPCA_DIR + "/{pca_panel}/smartpca.par",
        evec=PCA_SMARTPCA_DIR + "/{pca_panel}/smartpca.evec",
        eval=PCA_SMARTPCA_DIR + "/{pca_panel}/smartpca.eval"
    log:
        "logs/population_structure/pca/smartpca/{pca_panel}.smartpca.log"
    params:
        n_components=lambda wildcards: config["population_structure"]["pca"]["smartpca"]["n_components"],
        numoutlieriter=lambda wildcards: config["population_structure"]["pca"]["smartpca"]["numoutlieriter"],
        usenorm=lambda wildcards: config["population_structure"]["pca"]["smartpca"]["usenorm"],
        altnormstyle=lambda wildcards: config["population_structure"]["pca"]["smartpca"]["altnormstyle"],
        familynames=lambda wildcards: config["population_structure"]["pca"]["smartpca"]["familynames"],
        fastmode=lambda wildcards: config["population_structure"]["pca"]["smartpca"]["fastmode"],
        numchrom=lambda wildcards: config["population_structure"]["pca"]["smartpca"]["numchrom"]
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PCA_SMARTPCA_DIR}/{wildcards.pca_panel} logs/population_structure/pca/smartpca

        cat > {output.par} <<EOF
genotypename: {input.bed}
snpname: {input.bim}
indivname: {input.fam}
evecoutname: {output.evec}
evaloutname: {output.eval}
numoutevec: {params.n_components}
numoutlieriter: {params.numoutlieriter}
usenorm: {params.usenorm}
altnormstyle: {params.altnormstyle}
familynames: {params.familynames}
fastmode: {params.fastmode}
numchrom: {params.numchrom}
EOF

        export OMP_NUM_THREADS={threads}
        smartpca -p {output.par} > {log} 2>&1

        test -s {output.evec}
        test -s {output.eval}
        """


rule summarize_and_plot_joint_smartpca:
    input:
        evec=PCA_SMARTPCA_DIR + "/{pca_panel}/smartpca.evec",
        eval=PCA_SMARTPCA_DIR + "/{pca_panel}/smartpca.eval",
        annotations=PCA_SAMPLE_ANNOTATIONS,
        script="workflow/scripts/summarize_and_plot_smartpca.py"
    output:
        coordinates=PCA_SMARTPCA_DIR + "/{pca_panel}/pca_coordinates.tsv",
        eigenvalues=PCA_SMARTPCA_DIR + "/{pca_panel}/pca_eigenvalues.tsv",
        summary=PCA_SMARTPCA_DIR + "/{pca_panel}/pca_summary.tsv",
        pc12_png=PCA_SMARTPCA_DIR + "/{pca_panel}/pc1_pc2.png",
        pc12_pdf=PCA_SMARTPCA_DIR + "/{pca_panel}/pc1_pc2.pdf",
        pc23_png=PCA_SMARTPCA_DIR + "/{pca_panel}/pc2_pc3.png",
        pc23_pdf=PCA_SMARTPCA_DIR + "/{pca_panel}/pc2_pc3.pdf"
    params:
        panel=lambda wildcards: wildcards.pca_panel,
        expected_samples=lambda wildcards: config["population_structure"]["pca"]["expected_joint_samples"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --evec {input.evec} \
            --eval {input.eval} \
            --annotations {input.annotations} \
            --panel {params.panel} \
            --expected-samples {params.expected_samples} \
            --coordinates-out {output.coordinates} \
            --eigenvalues-out {output.eigenvalues} \
            --summary-out {output.summary} \
            --pc12-png {output.pc12_png} \
            --pc12-pdf {output.pc12_pdf} \
            --pc23-png {output.pc23_png} \
            --pc23-pdf {output.pc23_pdf}
        """


rule run_joint_smartpca_both_panels:
    input:
        PCA_SAMPLE_ANNOTATIONS,
        PCA_SAMPLE_COUNTS,
        expand(PCA_SMARTPCA_DIR + "/{pca_panel}/smartpca.par", pca_panel=["unmasked", "highld_masked"]),
        expand(PCA_SMARTPCA_DIR + "/{pca_panel}/smartpca.evec", pca_panel=["unmasked", "highld_masked"]),
        expand(PCA_SMARTPCA_DIR + "/{pca_panel}/smartpca.eval", pca_panel=["unmasked", "highld_masked"]),
        expand(PCA_SMARTPCA_DIR + "/{pca_panel}/pca_coordinates.tsv", pca_panel=["unmasked", "highld_masked"]),
        expand(PCA_SMARTPCA_DIR + "/{pca_panel}/pca_eigenvalues.tsv", pca_panel=["unmasked", "highld_masked"]),
        expand(PCA_SMARTPCA_DIR + "/{pca_panel}/pca_summary.tsv", pca_panel=["unmasked", "highld_masked"]),
        expand(PCA_SMARTPCA_DIR + "/{pca_panel}/pc1_pc2.png", pca_panel=["unmasked", "highld_masked"]),
        expand(PCA_SMARTPCA_DIR + "/{pca_panel}/pc1_pc2.pdf", pca_panel=["unmasked", "highld_masked"]),
        expand(PCA_SMARTPCA_DIR + "/{pca_panel}/pc2_pc3.png", pca_panel=["unmasked", "highld_masked"]),
        expand(PCA_SMARTPCA_DIR + "/{pca_panel}/pc2_pc3.pdf", pca_panel=["unmasked", "highld_masked"])
