# ============================================================
# Selection branch: production phased VCF inputs for saltiLASSI/LASSI
# ============================================================
#
# This module prepares and audits the exact marker panels that will be handed
# to lassip. It deliberately stops before computing HFS spectra or scan
# statistics.
#
# CT:
#   - converts the audited SHAPEIT2 HAPS output back to a phased VCF;
#   - restores the original source-VCF REF/ALT orientation;
#   - audits the exact post-SHAPEIT2 marker density against the pre-phasing
#     LD-derived 55.5-kb scale.
#
# CEU/TSI/IBS:
#   - subsets the published phased 1000G Phase 3 VCFs by population;
#   - applies the same population-specific SNP/MAF/missingness definition used
#     for the corresponding LD-calibration panel;
#   - requires exact variant-ID/order identity with that calibration BIM.

LASSI_INPUT_DIR = "results/selection/lassi/inputs"
LASSI_INPUT_AUDIT_DIR = "results/selection/lassi/input_audit"
LASSI_INPUT_CHROMS = [str(c) for c in range(1, 23)]
LASSI_REF_POPS = ["CEU", "TSI", "IBS"]


def lassi_ref_source_vcf(wildcards):
    return os.path.join(
        config["local_paths"]["onekg_phase3_vcf_dir"],
        config["datasets"]["reference_1kg_eur"]["phase3_vcf_filename_template"].format(
            chrom=wildcards.chrom
        ),
    )


def lassi_ref_source_index(wildcards):
    return lassi_ref_source_vcf(wildcards) + ".tbi"


rule selection_lassi_ct_popfile:
    input:
        samples="results/selection/phasing/audit/sample_order.txt"
    output:
        popfile=LASSI_INPUT_DIR + "/CT/CT.pop.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.popfile})
        awk 'NF {{print $1, "CT"}}' {input.samples} > {output.popfile}
        [ "$(wc -l < {output.popfile})" -eq 46 ]
        """


rule selection_lassi_reference_popfile:
    input:
        metadata="results/reference/1kg_eur/metadata/eur_samples.tsv"
    output:
        popfile=LASSI_INPUT_DIR + "/{population}/{population}.pop.txt",
        samples=LASSI_INPUT_DIR + "/{population}/{population}.samples.txt"
    params:
        expected=lambda wc: config["datasets"]["reference_1kg_eur"]["populations"]["expected_counts"][wc.population]
    wildcard_constraints:
        population="CEU|TSI|IBS"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.popfile})
        awk -F '\t' -v pop="{wildcards.population}" 'NR>1 && $2==pop {{print $1}}' {input.metadata} > {output.samples}
        awk -v pop="{wildcards.population}" 'NF {{print $1, pop}}' {output.samples} > {output.popfile}
        n=$(wc -l < {output.samples})
        [ "$n" -eq "{params.expected}" ] || {{ echo "ERROR: {wildcards.population}: expected {params.expected} samples, found $n" >&2; exit 1; }}
        """


rule selection_lassi_ct_phased_vcf:
    input:
        haps="results/selection/phasing/shapeit2/phased/chr{chrom}.ct.maf005.phased.haps.gz",
        sample="results/selection/phasing/shapeit2/phased/chr{chrom}.ct.maf005.phased.sample",
        source_vcf="results/selection/phasing/shapeit2/input/chr{chrom}.ct.maf005.vcf.gz",
        post_audit="results/selection/phasing/audit/postphasing_audit_summary.tsv",
        converter="workflow/scripts/shapeit_haps_to_vcf.py"
    output:
        vcf=LASSI_INPUT_DIR + "/CT/chr{chrom}.CT.maf005.phased.vcf.gz"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.converter}             --haps {input.haps}             --sample {input.sample}             --source-vcf {input.source_vcf}             --output {output.vcf}
        test -s {output.vcf}
        """


rule selection_lassi_reference_phased_vcf:
    input:
        vcf=lassi_ref_source_vcf,
        index=lassi_ref_source_index,
        samples=LASSI_INPUT_DIR + "/{population}/{population}.samples.txt"
    output:
        vcf=LASSI_INPUT_DIR + "/{population}/chr{chrom}.{population}.maf005.phased.vcf.gz"
    log:
        "logs/selection/lassi/input/{population}.chr{chrom}.log"
    wildcard_constraints:
        population="CEU|TSI|IBS"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf}) logs/selection/lassi/input

        bcftools view             -S {input.samples}             -m2 -M2 -v snps             -Ou {input.vcf} 2> {log}         | bcftools +fill-tags -Ou -- -t AF,F_MISSING 2>> {log}         | bcftools view             -i 'INFO/AF>=0.05 && INFO/AF<=0.95 && INFO/F_MISSING<=0.05'             -Oz -o {output.vcf} 2>> {log}

        test -s {output.vcf}
        """


rule audit_selection_lassi_ct_input:
    input:
        vcfs=expand(
            LASSI_INPUT_DIR + "/CT/chr{chrom}.CT.maf005.phased.vcf.gz",
            chrom=LASSI_INPUT_CHROMS,
        ),
        popfile=LASSI_INPUT_DIR + "/CT/CT.pop.txt",
        script="workflow/scripts/audit_lassi_input_vcf.py"
    output:
        summary=LASSI_INPUT_AUDIT_DIR + "/CT/input_summary.tsv",
        by_chrom=LASSI_INPUT_AUDIT_DIR + "/CT/input_by_chromosome.tsv",
        bim=LASSI_INPUT_AUDIT_DIR + "/CT/CT.phased.maf005.bim"
    params:
        expected=lambda wc: config["datasets"]["target"]["historical_post_qc_sample_count_to_verify"],
        maf=lambda wc: config["selection"]["lassi"]["marker_panel"]["maf_threshold"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script}             --population CT             --expected-samples {params.expected}             --pop-file {input.popfile}             --vcfs {input.vcfs}             --maf-threshold {params.maf}             --max-site-missing 0.05             --summary-out {output.summary}             --by-chrom-out {output.by_chrom}             --panel-bim-out {output.bim}
        """


rule audit_selection_lassi_reference_input:
    input:
        vcfs=lambda wc: [
            f"{LASSI_INPUT_DIR}/{wc.population}/chr{chrom}.{wc.population}.maf005.phased.vcf.gz"
            for chrom in LASSI_INPUT_CHROMS
        ],
        popfile=LASSI_INPUT_DIR + "/{population}/{population}.pop.txt",
        expected_bim="results/selection/ld_decay/reference/{population}/{population}.maf005.bim",
        script="workflow/scripts/audit_lassi_input_vcf.py"
    output:
        summary=LASSI_INPUT_AUDIT_DIR + "/{population}/input_summary.tsv",
        by_chrom=LASSI_INPUT_AUDIT_DIR + "/{population}/input_by_chromosome.tsv",
        bim=LASSI_INPUT_AUDIT_DIR + "/{population}/{population}.phased.maf005.bim"
    params:
        expected=lambda wc: config["datasets"]["reference_1kg_eur"]["populations"]["expected_counts"][wc.population],
        maf=lambda wc: config["selection"]["lassi"]["marker_panel"]["maf_threshold"]
    wildcard_constraints:
        population="CEU|TSI|IBS"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script}             --population {wildcards.population}             --expected-samples {params.expected}             --pop-file {input.popfile}             --vcfs {input.vcfs}             --maf-threshold {params.maf}             --max-site-missing 0.05             --require-maf             --expected-bim {input.expected_bim}             --summary-out {output.summary}             --by-chrom-out {output.by_chrom}             --panel-bim-out {output.bim}
        """


rule audit_selection_lassi_ct_phased_window_density:
    input:
        bim=LASSI_INPUT_AUDIT_DIR + "/CT/CT.phased.maf005.bim",
        anchor_bim="results/selection/ld_decay/panels/maf005/ct.maf005.bim",
        anchors="results/selection/ld_decay/anchors/maf005.anchors.txt",
        ld_summary="results/selection/ld_decay/summary/maf005.ld_decay_summary.tsv",
        script="workflow/scripts/derive_lassi_snp_window.py"
    output:
        anchor_counts=LASSI_INPUT_AUDIT_DIR + "/CT/phased_window_density/anchor_snp_counts.tsv",
        chrom_summary=LASSI_INPUT_AUDIT_DIR + "/CT/phased_window_density/chromosome_snp_window_summary.tsv",
        summary=LASSI_INPUT_AUDIT_DIR + "/CT/phased_window_density/lassi_window_summary.tsv"
    params:
        step_fraction=lambda wc: config["selection"]["lassi"]["windowing"]["salti_step_fraction"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.summary})
        primary=$(awk -F '\t' '$1=="first_crossing_bin_center_kb" {{print $2}}' {input.ld_summary})
        stable=$(awk -F '\t' '$1=="stable_crossing_bin_center_kb" {{print $2}}' {input.ld_summary})
        [ -n "$primary" ] && [ "$primary" != "NA" ] || {{ echo "ERROR: missing CT primary LD crossing" >&2; exit 1; }}
        [ -n "$stable" ] && [ "$stable" != "NA" ] || stable="$primary"

        python {input.script}             --bim {input.bim}             --anchor-bim {input.anchor_bim}             --anchors {input.anchors}             --primary-width-kb "$primary"             --stable-width-kb "$stable"             --step-fraction {params.step_fraction}             --anchor-counts-out {output.anchor_counts}             --chrom-summary-out {output.chrom_summary}             --summary-out {output.summary}
        """


rule audit_selection_lassi_input_panels:
    input:
        ct=LASSI_INPUT_AUDIT_DIR + "/CT/input_summary.tsv",
        ct_density=LASSI_INPUT_AUDIT_DIR + "/CT/phased_window_density/lassi_window_summary.tsv",
        refs=expand(
            LASSI_INPUT_AUDIT_DIR + "/{population}/input_summary.tsv",
            population=LASSI_REF_POPS,
        )
