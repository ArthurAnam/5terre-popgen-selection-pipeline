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
        audit="results/selection/lassi/preflight/crg100_resource_audit.tsv"
    output:
        ok="results/selection/lassi/preflight/crg100_resource_gate.ok"
    shell:
        r"""
        set -euo pipefail

        status=$(awk -F '\t' '$1=="overall_crg100_resource_status" {print $2; exit}' {input.audit})
        if [ "$status" != "PASS" ]; then
            echo "CRG100 resource gate: FAIL or missing overall PASS status" >&2
            awk -F '\t' 'NR>1 && $2=="FAIL" {print "  FAIL\t"$1"\t"$3}' {input.audit} >&2 || true
            exit 1
        fi

        printf 'PASS\n' > {output.ok}
        echo "CRG100 resource gate: PASS"
        """
