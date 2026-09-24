# ============================================================
# Pre-LASSI system / repository / local-file audit
# ============================================================

rule audit_pre_lassi_environment:
    input:
        phasing_audit="results/selection/phasing/shapeit2/audit/phasing_audit_summary.tsv",
        script="workflow/scripts/audit_pre_lassi_environment.py"
    output:
        summary="results/provenance/pre_lassi_environment_audit.tsv",
        details="results/provenance/pre_lassi_environment_details.txt",
        manifest="results/provenance/pre_lassi_critical_manifest.tsv",
        conda_list="results/provenance/conda_list.txt",
        conda_explicit="results/provenance/conda_explicit.txt",
        conda_history="results/provenance/conda_from_history.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script}
        test -s {output.summary}
        test -s {output.details}
        test -s {output.manifest}
        test -s {output.conda_list}
        test -s {output.conda_explicit}
        test -s {output.conda_history}
        """

# The diagnostic audit above always preserves a completed report. This separate
# gate is the blocking target to depend on before any production LASSI scan.
rule gate_pre_lassi_environment:
    input:
        summary="results/provenance/pre_lassi_environment_audit.tsv",
        script="workflow/scripts/gate_pre_lassi_environment.py"
    output:
        ok="results/provenance/pre_lassi_environment_gate.ok"
    shell:
        r"""
        set -euo pipefail
        python {input.script} {input.summary} {output.ok}
        test -s {output.ok}
        """
