# ============================================================
# Pinned lassip v1.2.1 software + method-configuration audit and blocking gate
# ============================================================

rule audit_selection_lassip_software:
    input:
        script="workflow/scripts/audit_lassip_software.py",
        config="config/config.yaml"
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
        audit="results/selection/lassi/preflight/lassip_software_audit.tsv"
    output:
        ok="results/selection/lassi/preflight/lassip_software_gate.ok"
    shell:
        r"""
        set -euo pipefail

        status=$(awk -F '\t' '$1=="overall_lassip_software_status" {print $2; exit}' {input.audit})
        if [ "$status" != "PASS" ]; then
            echo "lassip software gate: FAIL or missing overall PASS status" >&2
            awk -F '\t' 'NR>1 && $2=="FAIL" {print "  FAIL\t"$1"\t"$3}' {input.audit} >&2 || true
            exit 1
        fi

        printf 'PASS\n' > {output.ok}
        echo "lassip software gate: PASS"
        """
