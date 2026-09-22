# ============================================================
# Selection branch: convert physical LD decay to LASSI SNP window
# ============================================================

LASSI_WINDOW_DIR = "results/selection/lassi_window"

rule derive_lassi_snp_window:
    input:
        bim=SEL_LD_DIR + "/panels/maf005/ct.maf005.bim",
        anchors=SEL_LD_DIR + "/anchors/maf005.anchors.txt",
        script="workflow/scripts/derive_lassi_snp_window.py"
    output:
        anchor_counts=LASSI_WINDOW_DIR + "/anchor_snp_counts.tsv",
        chrom_summary=LASSI_WINDOW_DIR + "/chromosome_snp_window_summary.tsv",
        summary=LASSI_WINDOW_DIR + "/lassi_window_summary.tsv"
    params:
        primary_kb=lambda wc: config["selection"]["ld_decay"]["derived_parameters"]["selected_physical_decay_kb"],
        stable_kb=lambda wc: config["selection"]["ld_decay"]["derived_parameters"]["stable_diagnostic_kb"],
        step_fraction=lambda wc: config["selection"]["ld_decay"]["derived_parameters"]["shift_fraction_of_window"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {LASSI_WINDOW_DIR}
        python {input.script} \
            --bim {input.bim} \
            --anchors {input.anchors} \
            --primary-width-kb {params.primary_kb} \
            --stable-width-kb {params.stable_kb} \
            --step-fraction {params.step_fraction} \
            --anchor-counts-out {output.anchor_counts} \
            --chrom-summary-out {output.chrom_summary} \
            --summary-out {output.summary}
        """

rule audit_lassi_snp_window:
    input:
        LASSI_WINDOW_DIR + "/anchor_snp_counts.tsv",
        LASSI_WINDOW_DIR + "/chromosome_snp_window_summary.tsv",
        LASSI_WINDOW_DIR + "/lassi_window_summary.tsv"
