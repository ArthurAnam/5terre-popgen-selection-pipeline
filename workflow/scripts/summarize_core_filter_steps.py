#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import OrderedDict


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Count the sequential effect of the agreed core site filters: "
            "monomorphic sites, palindromic A/T or C/G SNPs, then site missingness >5%."
        )
    )
    p.add_argument("--out", required=True)
    return p.parse_args()


def num(text):
    if text in ("", ".", "NA", None):
        return None
    return float(text)


def palindromic(ref, alt):
    pair = {ref.upper(), alt.upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}


def main():
    args = parse_args()
    reader = csv.DictReader(sys.stdin, delimiter="\t")
    required = {"REF", "ALT", "MAF", "F_MISSING"}
    if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
        raise ValueError("Missing required input fields")

    counts = OrderedDict([
        ("starting_autosomal_sites", 0),
        ("removed_monomorphic", 0),
        ("after_monomorphic", 0),
        ("removed_palindromic_AT_CG", 0),
        ("after_palindromic", 0),
        ("removed_variant_missingness_gt_0.05", 0),
        ("after_variant_missingness", 0),
    ])

    after_mono = 0
    after_pal = 0
    after_miss = 0

    for row in reader:
        counts["starting_autosomal_sites"] += 1
        maf = num(row["MAF"])
        f_missing = num(row["F_MISSING"])

        if maf is not None and maf == 0:
            counts["removed_monomorphic"] += 1
            continue
        after_mono += 1

        if palindromic(row["REF"], row["ALT"]):
            counts["removed_palindromic_AT_CG"] += 1
            continue
        after_pal += 1

        if f_missing is not None and f_missing > 0.05:
            counts["removed_variant_missingness_gt_0.05"] += 1
            continue
        after_miss += 1

    counts["after_monomorphic"] = after_mono
    counts["after_palindromic"] = after_pal
    counts["after_variant_missingness"] = after_miss

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("step\tcount\n")
        for key, value in counts.items():
            handle.write(f"{key}\t{value}\n")


if __name__ == "__main__":
    main()
