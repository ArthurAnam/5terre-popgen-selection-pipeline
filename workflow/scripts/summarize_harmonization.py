#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def read_audit(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly one audit row, found {len(rows)}")
    row = rows[0]
    for key in (
        "CT",
        "EUR_QC",
        "same_position",
        "exact_match",
        "allele_mismatch",
        "CT_samples",
        "EUR_samples",
        "merged_samples",
    ):
        row[key] = int(row[key])
    return row


def count_pvar_variants(path):
    n = 0
    with open(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                n += 1
    return n


def count_psam_samples(path):
    n = 0
    with open(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def chrom_number(label):
    return int(label.removeprefix("chr"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audits", nargs="+", required=True)
    parser.add_argument("--pvar", required=True)
    parser.add_argument("--psam", required=True)
    parser.add_argument("--by-chromosome-out", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    rows = [read_audit(path) for path in args.audits]
    rows.sort(key=lambda row: chrom_number(row["chr"]))

    if [chrom_number(row["chr"]) for row in rows] != list(range(1, 23)):
        raise ValueError("Expected exactly chromosomes 1-22 in harmonization audits")

    fieldnames = [
        "chr",
        "CT",
        "EUR_QC",
        "same_position",
        "exact_match",
        "allele_mismatch",
        "CT_samples",
        "EUR_samples",
        "merged_samples",
    ]

    Path(args.by_chromosome_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.by_chromosome_out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        if row["same_position"] - row["exact_match"] != row["allele_mismatch"]:
            raise ValueError(
                f'{row["chr"]}: allele_mismatch is not same_position - exact_match'
            )

    sample_tuples = {
        (row["CT_samples"], row["EUR_samples"], row["merged_samples"]) for row in rows
    }
    if len(sample_tuples) != 1:
        raise ValueError("Sample counts are not identical across chromosomes")

    ct_samples, eur_samples, merged_samples = sample_tuples.pop()
    if merged_samples != ct_samples + eur_samples:
        raise ValueError("Merged sample count does not equal CT + EUR sample counts")

    totals = {
        key: sum(row[key] for row in rows)
        for key in ("CT", "EUR_QC", "same_position", "exact_match", "allele_mismatch")
    }

    pvar_variants = count_pvar_variants(args.pvar)
    psam_samples = count_psam_samples(args.psam)

    if pvar_variants != totals["exact_match"]:
        raise ValueError(
            f"PVAR has {pvar_variants} variants but audits contain "
            f'{totals["exact_match"]} exact matches'
        )
    if psam_samples != merged_samples:
        raise ValueError(
            f"PSAM has {psam_samples} samples but chromosome audits contain "
            f"{merged_samples}"
        )

    with open(args.summary_out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["target_qc_variants", totals["CT"]])
        writer.writerow(["eur_qc_nonpalindromic_variants", totals["EUR_QC"]])
        writer.writerow(["shared_positions", totals["same_position"]])
        writer.writerow(["exact_chr_pos_ref_alt_matches", totals["exact_match"]])
        writer.writerow(["same_position_allele_mismatches", totals["allele_mismatch"]])
        writer.writerow(["target_samples", ct_samples])
        writer.writerow(["eur_samples", eur_samples])
        writer.writerow(["harmonized_samples", merged_samples])
        writer.writerow(["harmonized_pgen_variants", pvar_variants])
        writer.writerow(["match_policy", "exact_CHR:POS:REF:ALT_only"])
        writer.writerow(["allele_rescue_or_strand_flip", "NO"])


if __name__ == "__main__":
    main()
