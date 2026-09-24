# ============================================================
# Pinned lassip v1.2.1 software audit and blocking gate
# ============================================================

rule audit_selection_lassip_software:
    input:
        script="workflow/scripts/audit_lassip_software.py"
    output:
        audit="results/selection/lassi/preflight/lassip_software_audit.tsv"
    shell:
        r"""
        set -euo pipefail
        python {input.script}
        test -s {output.audit}
        """

rule gate_selection_lassip_software:
    input:
        audit="results/selection/lassi/preflight/lassip_software_audit.tsv",
        script="workflow/scripts/gate_lassip_software.py"
    output:
        ok="results/selection/lassi/preflight/lassip_software_gate.ok"
    shell:
        r"""
        set -euo pipefail
        python {input.script} {input.audit} {output.ok}
        test -s {output.ok}
        """
