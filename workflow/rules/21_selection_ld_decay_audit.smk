# ============================================================
# Selection branch: CT LD-decay audit for LASSI window calibration
# ============================================================

SEL_LD_DIR = "results/selection/ld_decay"
SEL_LD_TAGS = ["maf005", "maf001"]
SEL_LD_MAF = {"maf005": 0.05, "maf001": 0.01}
SEL_LD_LABEL = {"maf005": "MAF>=0.05", "maf001": "MAF>=0.01"}
SEL_LD_VCF = "results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz"


rule selection_ld_make_panel:
    input:
        vcf=SEL_LD_VCF,
        tbi=SEL_LD_VCF + ".tbi"
    output:
        bed=SEL_LD_DIR + "/panels/{maf_tag}/ct.{maf_tag}.bed",
        bim=SEL_LD_DIR + "/panels/{maf_tag}/ct.{maf_tag}.bim",
        fam=SEL_LD_DIR + "/panels/{maf_tag}/ct.{maf_tag}.fam"
    params:
        maf=lambda wc: SEL_LD_MAF[wc.maf_tag],
        prefix=lambda wc: f"{SEL_LD_DIR}/panels/{wc.maf_tag}/ct.{wc.maf_tag}"
    log:
        "logs/selection/ld_decay/{maf_tag}.panel.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bed}) logs/selection/ld_decay

        plink2 \
            --threads {threads} \
            --vcf {input.vcf} \
            --maf {params.maf} \
            --set-all-var-ids '@:#:$r:$a' \
            --make-bed \
            --out {params.prefix} \
            > {log} 2>&1
        """


rule selection_ld_make_anchors:
    input:
        bim=SEL_LD_DIR + "/panels/{maf_tag}/ct.{maf_tag}.bim",
        script="workflow/scripts/make_selection_ld_anchors.py"
    output:
        anchors=SEL_LD_DIR + "/anchors/{maf_tag}.anchors.txt",
        summary=SEL_LD_DIR + "/anchors/{maf_tag}.anchor_summary.tsv"
    params:
        n=lambda wc: config["selection"]["ld_decay"]["anchor_sampling"]["anchors_per_chromosome"],
        edge=lambda wc: config["selection"]["ld_decay"]["anchor_sampling"]["edge_margin_kb"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --bim {input.bim} \
            --anchors-per-chromosome {params.n} \
            --edge-margin-kb {params.edge} \
            --anchors-out {output.anchors} \
            --summary-out {output.summary}
        """


rule selection_ld_compute_pairs:
    input:
        bed=SEL_LD_DIR + "/panels/{maf_tag}/ct.{maf_tag}.bed",
        bim=SEL_LD_DIR + "/panels/{maf_tag}/ct.{maf_tag}.bim",
        fam=SEL_LD_DIR + "/panels/{maf_tag}/ct.{maf_tag}.fam",
        anchors=SEL_LD_DIR + "/anchors/{maf_tag}.anchors.txt"
    output:
        ld=SEL_LD_DIR + "/pairs/{maf_tag}.ld.gz"
    params:
        prefix=lambda wc: f"{SEL_LD_DIR}/pairs/{wc.maf_tag}",
        max_kb=lambda wc: config["selection"]["ld_decay"]["anchor_sampling"]["max_pair_distance_kb"]
    log:
        "logs/selection/ld_decay/{maf_tag}.r2.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SEL_LD_DIR}/pairs logs/selection/ld_decay

        plink \
            --bfile {SEL_LD_DIR}/panels/{wildcards.maf_tag}/ct.{wildcards.maf_tag} \
            --r2 gz \
            --ld-snp-list {input.anchors} \
            --ld-window 99999999 \
            --ld-window-kb {params.max_kb} \
            --ld-window-r2 0 \
            --out {params.prefix} \
            > {log} 2>&1
        """


rule selection_ld_summarize:
    input:
        ld=SEL_LD_DIR + "/pairs/{maf_tag}.ld.gz",
        script="workflow/scripts/summarize_selection_ld_decay.py"
    output:
        bins=SEL_LD_DIR + "/summary/{maf_tag}.ld_decay_bins.tsv",
        summary=SEL_LD_DIR + "/summary/{maf_tag}.ld_decay_summary.tsv",
        png=SEL_LD_DIR + "/summary/{maf_tag}.ld_decay.png",
        pdf=SEL_LD_DIR + "/summary/{maf_tag}.ld_decay.pdf",
        readme=SEL_LD_DIR + "/summary/{maf_tag}.README.txt"
    params:
        label=lambda wc: SEL_LD_LABEL[wc.maf_tag],
        baseline=lambda wc: config["selection"]["ld_decay"]["baseline_distance_kb"],
        half=lambda wc: config["selection"]["ld_decay"]["baseline_half_width_kb"],
        fraction=lambda wc: config["selection"]["ld_decay"]["decay_fraction"],
        bin_width=lambda wc: config["selection"]["ld_decay"]["bin_width_kb"],
        max_kb=lambda wc: config["selection"]["ld_decay"]["anchor_sampling"]["max_pair_distance_kb"],
        stable=lambda wc: config["selection"]["ld_decay"]["stability_diagnostic_consecutive_bins"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SEL_LD_DIR}/summary

        python {input.script} \
            --ld-gz {input.ld} \
            --maf-label '{params.label}' \
            --baseline-kb {params.baseline} \
            --baseline-half-width-kb {params.half} \
            --decay-fraction {params.fraction} \
            --bin-width-kb {params.bin_width} \
            --max-distance-kb {params.max_kb} \
            --stable-bins {params.stable} \
            --bins-out {output.bins} \
            --summary-out {output.summary} \
            --plot-png {output.png} \
            --plot-pdf {output.pdf} \
            --readme-out {output.readme}
        """


rule audit_selection_ld_decay:
    input:
        SEL_LD_DIR + "/anchors/maf005.anchor_summary.tsv",
        SEL_LD_DIR + "/summary/maf005.ld_decay_bins.tsv",
        SEL_LD_DIR + "/summary/maf005.ld_decay_summary.tsv",
        SEL_LD_DIR + "/summary/maf005.ld_decay.png",
        SEL_LD_DIR + "/summary/maf005.ld_decay.pdf",
        SEL_LD_DIR + "/summary/maf005.README.txt"


rule audit_selection_ld_decay_sensitivity:
    input:
        SEL_LD_DIR + "/anchors/maf001.anchor_summary.tsv",
        SEL_LD_DIR + "/summary/maf001.ld_decay_bins.tsv",
        SEL_LD_DIR + "/summary/maf001.ld_decay_summary.tsv",
        SEL_LD_DIR + "/summary/maf001.ld_decay.png",
        SEL_LD_DIR + "/summary/maf001.ld_decay.pdf",
        SEL_LD_DIR + "/summary/maf001.README.txt"
