#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def as_int(row, key):
    return int(row[key])

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--audits", nargs="+", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--chrom-out", required=True)
    args = p.parse_args()

    rows = []
    for path in args.audits:
        with open(path, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows.extend(reader)

    rows.sort(key=lambda r: int(r["chromosome"]))
    if len(rows) != 22:
        raise SystemExit(f"ERROR: expected 22 chromosome audit rows, found {len(rows)}")

    Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)

    fields = list(rows[0].keys())
    with open(args.chrom_out, "w", encoding="utf-8", newline="") as out:
        w = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    ct_bcf = sum(as_int(r, "ct_maf005_bcftools") for r in rows)
    ct_plink = sum(as_int(r, "ct_maf005_plink") for r in rows)
    eur_exact = sum(as_int(r, "exact_eur_polymorphic") for r in rows)
    eur_same = sum(as_int(r, "same_position_eur_polymorphic") for r in rows)
    eur_ref = sum(as_int(r, "eur_polymorphic_biallelic_snps") for r in rows)
    mismatches = [r for r in rows if r["ct_panel_count_match"] != "PASS"]

    with open(args.summary_out, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["autosomes", len(rows)])
        w.writerow(["ct_maf005_bcftools_total", ct_bcf])
        w.writerow(["ct_maf005_plink_total", ct_plink])
        w.writerow(["ct_panel_count_match_all_chromosomes", "PASS" if not mismatches else "FAIL"])
        w.writerow(["eur_polymorphic_biallelic_snps_total", eur_ref])
        w.writerow(["exact_eur_polymorphic_total", eur_exact])
        w.writerow(["fraction_ct_maf005_exact_eur_polymorphic", eur_exact / ct_bcf if ct_bcf else 0.0])
        w.writerow(["same_position_eur_polymorphic_total", eur_same])
        w.writerow(["same_position_allele_mismatch_eur_polymorphic", eur_same - eur_exact])
        w.writerow(["reference_scope", "1000G_EUR_FIXED"])
        w.writerow(["effective_population_size", 11418])
        w.writerow(["production_phasing_started", "NO"])

if __name__ == "__main__":
    main()
