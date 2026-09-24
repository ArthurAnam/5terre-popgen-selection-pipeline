#!/usr/bin/env python3
import csv
import sys
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit("usage: gate_crg100_resource.py AUDIT_TSV OK_SENTINEL")

audit = Path(sys.argv[1])
ok = Path(sys.argv[2])

with audit.open(newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

overall = next((r for r in rows if r["check"] == "overall_crg100_resource_status"), None)
if overall is None:
    print("CRG100 resource gate: missing overall status", file=sys.stderr)
    raise SystemExit(2)

if overall["status"] != "PASS":
    print(f"CRG100 resource gate: FAIL ({overall['value']})", file=sys.stderr)
    for r in rows:
        if r["status"] == "FAIL":
            print(f"  FAIL\t{r['check']}\t{r['value']}", file=sys.stderr)
    raise SystemExit(1)

ok.parent.mkdir(parents=True, exist_ok=True)
ok.write_text("PASS\n")
print(f"CRG100 resource gate: PASS; wrote {ok}")
