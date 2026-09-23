#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path

def nonempty_lines(path):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip())

def chrom_from_name(path):
    name = Path(path).name
    return name.split("chr", 1)[1].split(".", 1)[0]

def read_text(path):
    return Path(path).read_text(encoding="utf-8", errors="replace")

def extract_int(text, pattern, label, path):
    m = re.search(pattern, text, flags=re.MULTILINE)
    if not m:
        raise SystemExit(f"ERROR: could not parse {label} from {path}")
    return int(m.group(1))

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--preflight", required=True)
    p.add_argument("--exclusions", nargs="+", required=True)
    p.add_argument("--initial-logs", nargs="+", required=True)
    p.add_argument("--postcheck-logs", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--chrom-out", required=True)
    args = p.parse_args()

    with open(args.preflight, "r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.sort(key=lambda r: int(r["chromosome"]))

    groups = [args.exclusions, args.initial_logs, args.postcheck_logs]
    if len(rows) != 22 or any(len(g) != 22 for g in groups):
        raise SystemExit("ERROR: expected 22 chromosomes for preflight, exclusions, and SHAPEIT2 logs")

    exc_lines = {chrom_from_name(p): nonempty_lines(p) for p in args.exclusions}
    initial = {chrom_from_name(p): p for p in args.initial_logs}
    post = {chrom_from_name(p): p for p in args.postcheck_logs}

    out_rows = []
    for row in rows:
        chrom = row["chromosome"]
        n_input = int(row["ct_maf005_bcftools"])

        itxt = read_text(initial[chrom])
        ptxt = read_text(post[chrom])

        n_missing = extract_int(
            itxt, r"#Missing sites in reference panel\s*=\s*(\d+)",
            "missing sites", initial[chrom]
        )
        n_misaligned = extract_int(
            itxt, r"#Misaligned sites between panels\s*=\s*(\d+)",
            "misaligned sites", initial[chrom]
        )

        n_retained = extract_int(
            ptxt, r"\*\s+(\d+) SNPs included\s*$",
            "post-exclusion included SNPs", post[chrom]
        )
        n_excluded = extract_int(
            ptxt, r"\*\s+(\d+) SNPs excluded\s*$",
            "post-exclusion excluded SNPs", post[chrom]
        )

        if n_retained + n_excluded != n_input:
            raise SystemExit(
                f"ERROR: chr{chrom}: retained + excluded != input "
                f"({n_retained}+{n_excluded}!={n_input})"
            )
        if n_missing + n_misaligned != n_excluded:
            raise SystemExit(
                f"ERROR: chr{chrom}: initial missing + misaligned != actual excluded "
                f"({n_missing}+{n_misaligned}!={n_excluded})"
            )

        out_rows.append({
            "chromosome": chrom,
            "ct_maf005_snps": n_input,
            "missing_in_eur_reference": n_missing,
            "misaligned_to_eur_reference": n_misaligned,
            "shapeit2_actual_excluded": n_excluded,
            "exclude_file_nonempty_lines": exc_lines[chrom],
            "retained_for_phasing": n_retained,
            "retained_fraction": n_retained / n_input if n_input else 0.0,
            "post_exclusion_check": "PASS",
        })

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.chrom_out, "w", encoding="utf-8", newline="") as out:
        fields = [
            "chromosome", "ct_maf005_snps", "missing_in_eur_reference",
            "misaligned_to_eur_reference", "shapeit2_actual_excluded",
            "exclude_file_nonempty_lines", "retained_for_phasing",
            "retained_fraction", "post_exclusion_check"
        ]
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    total_input = sum(r["ct_maf005_snps"] for r in out_rows)
    total_missing = sum(r["missing_in_eur_reference"] for r in out_rows)
    total_misaligned = sum(r["misaligned_to_eur_reference"] for r in out_rows)
    total_exc = sum(r["shapeit2_actual_excluded"] for r in out_rows)
    total_exc_lines = sum(r["exclude_file_nonempty_lines"] for r in out_rows)
    total_ret = sum(r["retained_for_phasing"] for r in out_rows)

    with open(args.out, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["metric","value"])
        w.writerow(["autosomes",22])
        w.writerow(["ct_maf005_snps",total_input])
        w.writerow(["missing_in_eur_reference",total_missing])
        w.writerow(["misaligned_to_eur_reference",total_misaligned])
        w.writerow(["shapeit2_actual_excluded",total_exc])
        w.writerow(["exclude_file_nonempty_lines",total_exc_lines])
        w.writerow(["retained_for_phasing",total_ret])
        w.writerow(["retained_fraction",total_ret/total_input if total_input else 0.0])
        w.writerow(["minimum_chromosome_retained_fraction",min(r["retained_fraction"] for r in out_rows)])
        w.writerow(["maximum_chromosome_retained_fraction",max(r["retained_fraction"] for r in out_rows)])
        w.writerow(["post_exclusion_shapeit2_check_all_chromosomes","PASS"])
        w.writerow(["review_status","REVIEW_BEFORE_PRODUCTION_PHASING"])

if __name__ == "__main__":
    main()
