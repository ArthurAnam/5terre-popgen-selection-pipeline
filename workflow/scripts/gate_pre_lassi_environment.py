#!/usr/bin/env python3
import csv
import sys
from pathlib import Path


def die(message, code=2):
    print(message, file=sys.stderr)
    sys.exit(code)


if len(sys.argv) != 3:
    die("usage: gate_pre_lassi_environment.py AUDIT_TSV OK_SENTINEL")

summary = Path(sys.argv[1])
ok = Path(sys.argv[2])

if not summary.is_file() or summary.stat().st_size == 0:
    die(f"pre-LASSI gate: missing/empty audit summary: {summary}")

with summary.open(newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

if not rows or not {"check", "status", "value"}.issubset(rows[0]):
    die(f"pre-LASSI gate: malformed audit summary: {summary}")

overall = next(
    (row for row in rows if row["check"] == "overall_pre_lassi_environment_status"),
    None,
)
if overall is None:
    die("pre-LASSI gate: overall_pre_lassi_environment_status is missing")

failures = [row for row in rows if row["status"] == "FAIL"]

if overall["status"] != "PASS":
    print(
        f"pre-LASSI gate: FAIL ({overall['value']}); diagnostic report retained at {summary}",
        file=sys.stderr,
    )
    for row in failures:
        print(f"  FAIL\t{row['check']}\t{row['value']}", file=sys.stderr)
    sys.exit(1)

if failures:
    die(
        "pre-LASSI gate: inconsistent audit: overall status is PASS but FAIL rows are present"
    )

ok.parent.mkdir(parents=True, exist_ok=True)
ok.write_text("PASS\n")
print(f"pre-LASSI gate: PASS; wrote {ok}")
