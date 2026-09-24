# ============================================================
# Selection branch: population-specific LD decay for all 1000G EUR
# and LASSI-window calibration for CEU, TSI and IBS
# ============================================================

SELECTION_REF_POPS = ["CEU", "FIN", "GBR", "IBS", "TSI"]
SELECTION_LASSI_REF_POPS = ["CEU", "TSI", "IBS"]
SELECTION_REF_LD_DIR = "results/selection/ld_decay/reference"
SELECTION_REF_WIN_DIR = "results/selection/lassi_window/reference"


def selection_ref_vcfs(wildcards):
    return [
        os.path.join(
            config["local_paths"]["onekg_phase3_vcf_dir"],
            config["datasets"]["reference_1kg_eur"]["phase3_vcf_filename_template"].format(chrom=chrom),
        )
        for chrom in HARMONIZATION_CHROMS
    ]


def selection_ref_vcf_indexes(wildcards):
    return [path + ".tbi" for path in selection_ref_vcfs(wildcards)]


rule selection_reference_sample_list:
    input:
        metadata="results/reference/1kg_eur/metadata/eur_samples.tsv"
    output:
        samples=SELECTION_REF_LD_DIR + "/{population}/samples.txt"
    params:
        expected=lambda wc: config["datasets"]["reference_1kg_eur"]["populations"]["expected_counts"][wc.population]
    wildcard_constraints:
        population="CEU|FIN|GBR|IBS|TSI"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.samples})
        awk -F '\t' -v pop="{wildcards.population}" 'NR>1 && $2==pop {{print $1}}' {input.metadata} > {output.samples}
        n=$(wc -l < {output.samples})
        [ "$n" -eq "{params.expected}" ] || {{ echo "ERROR: {wildcards.population} sample count $n != expected {params.expected}" >&2; exit 1; }}
        """


rule selection_reference_maf005_panel:
    input:
        vcfs=selection_ref_vcfs,
        indexes=selection_ref_vcf_indexes,
        samples=SELECTION_REF_LD_DIR + "/{population}/samples.txt"
    output:
        bed=SELECTION_REF_LD_DIR + "/{population}/{population}.maf005.bed",
        bim=SELECTION_REF_LD_DIR + "/{population}/{population}.maf005.bim",
        fam=SELECTION_REF_LD_DIR + "/{population}/{population}.maf005.fam"
    params:
        prefix=lambda wc: f"{SELECTION_REF_LD_DIR}/{wc.population}/{wc.population}.maf005"
    log:
        "logs/selection/ld_decay/reference/{population}.panel.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bed}) logs/selection/ld_decay/reference
        tmpdir="{SELECTION_REF_LD_DIR}/{wildcards.population}/.tmp_panel.$$"
        mkdir -p "$tmpdir"
        trap 'rm -rf "$tmpdir"' EXIT

        i=0
        for vcf in {input.vcfs}; do
            i=$((i+1))
            out="$tmpdir/chr$i"
            filtered="$tmpdir/chr$i.nonpal.maf005.vcf.gz"

            bcftools view \
                -S {input.samples} \
                -m2 -M2 -v snps \
                -Ou "$vcf" 2>> {log} \
            | bcftools +fill-tags -Ou -- -t AF,F_MISSING 2>> {log} \
            | bcftools view \
                -e 'INFO/F_MISSING>0.05 || INFO/AF<0.05 || INFO/AF>0.95 || ((REF="A" && ALT="T") || (REF="T" && ALT="A") || (REF="C" && ALT="G") || (REF="G" && ALT="C"))' \
                -Oz -o "$filtered" 2>> {log}

            plink2 \
                --threads {threads} \
                --vcf "$filtered" \
                --snps-only just-acgt \
                --max-alleles 2 \
                --set-all-var-ids '@:#:$r:$a' \
                --make-bed \
                --out "$out" \
                >> {log} 2>&1
        done

        first="$tmpdir/chr1"
        : > "$tmpdir/merge_list.txt"
        for i in $(seq 2 22); do
            printf '%s\n' "$tmpdir/chr$i" >> "$tmpdir/merge_list.txt"
        done

        plink \
            --bfile "$first" \
            --merge-list "$tmpdir/merge_list.txt" \
            --make-bed \
            --out {params.prefix} \
            >> {log} 2>&1
        """


rule selection_reference_ld_anchors:
    input:
        bim=SELECTION_REF_LD_DIR + "/{population}/{population}.maf005.bim",
        script="workflow/scripts/make_selection_ld_anchors.py"
    output:
        anchors=SELECTION_REF_LD_DIR + "/{population}/{population}.anchors.txt",
        summary=SELECTION_REF_LD_DIR + "/{population}/{population}.anchor_summary.tsv"
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


rule selection_reference_ld_pairs:
    input:
        bed=SELECTION_REF_LD_DIR + "/{population}/{population}.maf005.bed",
        bim=SELECTION_REF_LD_DIR + "/{population}/{population}.maf005.bim",
        fam=SELECTION_REF_LD_DIR + "/{population}/{population}.maf005.fam",
        anchors=SELECTION_REF_LD_DIR + "/{population}/{population}.anchors.txt"
    output:
        ld=SELECTION_REF_LD_DIR + "/{population}/{population}.ld.gz"
    params:
        prefix=lambda wc: f"{SELECTION_REF_LD_DIR}/{wc.population}/{wc.population}",
        max_kb=lambda wc: config["selection"]["ld_decay"]["anchor_sampling"]["max_pair_distance_kb"]
    log:
        "logs/selection/ld_decay/reference/{population}.r2.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        plink \
            --bfile {SELECTION_REF_LD_DIR}/{wildcards.population}/{wildcards.population}.maf005 \
            --r2 gz \
            --ld-snp-list {input.anchors} \
            --ld-window 99999999 \
            --ld-window-kb {params.max_kb} \
            --ld-window-r2 0 \
            --out {params.prefix} \
            > {log} 2>&1
        """


rule selection_reference_ld_summary:
    input:
        ld=SELECTION_REF_LD_DIR + "/{population}/{population}.ld.gz",
        script="workflow/scripts/summarize_selection_ld_decay.py"
    output:
        bins=SELECTION_REF_LD_DIR + "/{population}/{population}.ld_decay_bins.tsv",
        summary=SELECTION_REF_LD_DIR + "/{population}/{population}.ld_decay_summary.tsv",
        png=SELECTION_REF_LD_DIR + "/{population}/{population}.ld_decay.png",
        pdf=SELECTION_REF_LD_DIR + "/{population}/{population}.ld_decay.pdf",
        readme=SELECTION_REF_LD_DIR + "/{population}/{population}.README.txt"
    params:
        label=lambda wc: f"{wc.population} MAF>=0.05",
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


rule selection_reference_lassi_window:
    input:
        bim=SELECTION_REF_LD_DIR + "/{population}/{population}.maf005.bim",
        anchors=SELECTION_REF_LD_DIR + "/{population}/{population}.anchors.txt",
        ld_summary=SELECTION_REF_LD_DIR + "/{population}/{population}.ld_decay_summary.tsv",
        script="workflow/scripts/derive_lassi_snp_window.py"
    output:
        anchor_counts=SELECTION_REF_WIN_DIR + "/{population}/anchor_snp_counts.tsv",
        chrom_summary=SELECTION_REF_WIN_DIR + "/{population}/chromosome_snp_window_summary.tsv",
        summary=SELECTION_REF_WIN_DIR + "/{population}/lassi_window_summary.tsv"
    wildcard_constraints:
        population="CEU|TSI|IBS"
    params:
        step_fraction=lambda wc: config["selection"]["ld_decay"]["derived_parameters"]["shift_fraction_of_window"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SELECTION_REF_WIN_DIR}/{wildcards.population}
        primary=$(awk -F '\t' '$1=="first_crossing_bin_center_kb" {{print $2}}' {input.ld_summary})
        stable=$(awk -F '\t' '$1=="stable_crossing_bin_center_kb" {{print $2}}' {input.ld_summary})
        [ -n "$primary" ] && [ "$primary" != "NA" ] || {{ echo "ERROR: missing primary LD crossing" >&2; exit 1; }}
        [ -n "$stable" ] && [ "$stable" != "NA" ] || stable="$primary"

        python {input.script} \
            --bim {input.bim} \
            --anchors {input.anchors} \
            --primary-width-kb "$primary" \
            --stable-width-kb "$stable" \
            --step-fraction {params.step_fraction} \
            --anchor-counts-out {output.anchor_counts} \
            --chrom-summary-out {output.chrom_summary} \
            --summary-out {output.summary}
        """


rule audit_selection_reference_ld_windows:
    input:
        expand(SELECTION_REF_LD_DIR + "/{population}/{population}.ld_decay_summary.tsv", population=SELECTION_REF_POPS),
        expand(SELECTION_REF_LD_DIR + "/{population}/{population}.ld_decay.png", population=SELECTION_REF_POPS),
        expand(SELECTION_REF_LD_DIR + "/{population}/{population}.ld_decay.pdf", population=SELECTION_REF_POPS),
        expand(SELECTION_REF_WIN_DIR + "/{population}/lassi_window_summary.tsv", population=SELECTION_LASSI_REF_POPS),
        expand(SELECTION_REF_WIN_DIR + "/{population}/chromosome_snp_window_summary.tsv", population=SELECTION_LASSI_REF_POPS)
