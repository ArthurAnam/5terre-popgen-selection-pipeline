#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def nonempty_lines(path):
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                n += 1
    return n

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--preflight", required=True)
    p.add_argument("--exclusions", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--chrom-out", required=True)
    args = p.parse_args()

    with open(args.preflight, "r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.sort(key=lambda r: int(r["chromosome"]))
    if len(rows) != 22 or len(args.exclusions) != 22:
        raise SystemExit("ERROR: expected 22 preflight rows and 22 SHAPEIT2 exclusion files")

    exc_by_chr = {}
    for path in args.exclusions:
        name = Path(path).name
        chrom = name.split("chr", 1)[1].split(".", 1)[0]
        exc_by_chr[chrom] = nonempty_lines(path)

    out_rows = []
    for row in rows:
        chrom = row["chromosome"]
        n_input = int(row["ct_maf005_bcftools"])
        n_exc = exc_by_chr[chrom]
        retained = n_input - n_exc
        if retained < 0:
            raise SystemExit(f"ERROR: exclusions exceed input count on chr{chrom}")
        out_rows.append({
            "chromosome": chrom,
            "ct_maf005_snps": n_input,
            "shapeit2_check_exclusions": n_exc,
            "retained_for_phasing": retained,
            "retained_fraction": retained / n_input if n_input else 0.0,
        })

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.chrom_out, "w", encoding="utf-8", newline="") as out:
        fields = ["chromosome","ct_maf005_snps","shapeit2_check_exclusions","retained_for_phasing","retained_fraction"]
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    total_input = sum(r["ct_maf005_snps"] for r in out_rows)
    total_exc = sum(r["shapeit2_check_exclusions"] for r in out_rows)
    total_ret = sum(r["retained_for_phasing"] for r in out_rows)
    with open(args.out, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["metric","value"])
        w.writerow(["autosomes",22])
        w.writerow(["ct_maf005_snps",total_input])
        w.writerow(["shapeit2_check_exclusions",total_exc])
        w.writerow(["retained_for_phasing",total_ret])
        w.writerow(["retained_fraction",total_ret/total_input if total_input else 0.0])
        w.writerow(["minimum_chromosome_retained_fraction",min(r["retained_fraction"] for r in out_rows)])
        w.writerow(["maximum_chromosome_retained_fraction",max(r["retained_fraction"] for r in out_rows)])
        w.writerow(["review_status","REVIEW_BEFORE_PRODUCTION_PHASING"])

if __name__ == "__main__":
    main()
