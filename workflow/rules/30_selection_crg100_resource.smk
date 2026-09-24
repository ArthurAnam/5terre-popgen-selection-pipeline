# ============================================================
# GRCh37/hg19 CRG100 mappability resource audit and blocking gate
# ============================================================

rule audit_selection_crg100_resource:
    input:
        script="workflow/scripts/audit_crg100_resource.py",
        config="config/config.yaml"
    output:
        audit="results/selection/lassi/preflight/crg100_resource_audit.tsv"
    shell:
        r"""
        set -euo pipefail
        python {input.script}
        test -s {output.audit}
        """

rule gate_selection_crg100_resource:
    input:
        audit="results/selection/lassi/preflight/crg100_resource_audit.tsv",
        script="workflow/scripts/gate_crg100_resource.py"
    output:
        ok="results/selection/lassi/preflight/crg100_resource_gate.ok"
    shell:
        r"""
        set -euo pipefail
        python {input.script} {input.audit} {output.ok}
        test -s {output.ok}
        """
