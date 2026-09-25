#!/usr/bin/env python3

import argparse
from pathlib import Path


def count_nonempty(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def count_pvar(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def count_psam(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-pvar", required=True)
    parser.add_argument("--input-psam", required=True)
    parser.add_argument("--maf-snplist", required=True)
    parser.add_argument("--prune-in", required=True)
    parser.add_argument("--prune-out", required=True)
    parser.add_argument("--final-pvar", required=True)
    parser.add_argument("--final-psam", required=True)
    parser.add_argument("--maf-threshold", type=float, required=True)
    parser.add_argument("--maf-scope", required=True)
    parser.add_argument("--window-snps", type=int, required=True)
    parser.add_argument("--step-snps", type=int, required=True)
    parser.add_argument("--r2-threshold", type=float, required=True)
    parser.add_argument("--indep-order", type=int, required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    parser.add_argument("--expected-input-variants", type=int, required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    input_variants = count_pvar(Path(args.input_pvar))
    input_samples = count_psam(Path(args.input_psam))
    maf_variants = count_nonempty(Path(args.maf_snplist))
    prune_in = count_nonempty(Path(args.prune_in))
    prune_out = count_nonempty(Path(args.prune_out))
    final_variants = count_pvar(Path(args.final_pvar))
    final_samples = count_psam(Path(args.final_psam))

    if input_variants != args.expected_input_variants:
        raise SystemExit(
            f"ERROR: expected {args.expected_input_variants} harmonized variants, "
            f"found {input_variants}"
        )
    if input_samples != args.expected_samples:
        raise SystemExit(
            f"ERROR: expected {args.expected_samples} harmonized samples, "
            f"found {input_samples}"
        )
    if prune_in + prune_out != maf_variants:
        raise SystemExit(
            "ERROR: prune.in + prune.out does not equal the MAF-passing marker count "
            f"({prune_in} + {prune_out} != {maf_variants})"
        )
    if final_variants != prune_in:
        raise SystemExit(
            "ERROR: final pruned PVAR count does not equal prune.in count "
            f"({final_variants} != {prune_in})"
        )
    if final_samples != input_samples:
        raise SystemExit(
            "ERROR: PCA preprocessing unexpectedly changed the sample count "
            f"({final_samples} != {input_samples})"
        )

    removed_maf = input_variants - maf_variants
    removed_ld = prune_out

    rows = [
        ("harmonized_input_variants", input_variants),
        ("harmonized_input_samples", input_samples),
        ("maf_threshold", args.maf_threshold),
        ("maf_scope", args.maf_scope),
        ("variants_after_maf", maf_variants),
        ("variants_removed_by_maf", removed_maf),
        ("ld_window_snps", args.window_snps),
        ("ld_step_snps", args.step_snps),
        ("ld_r2_threshold", args.r2_threshold),
        ("plink2_indep_order", args.indep_order),
        ("variants_removed_by_ld_pruning", removed_ld),
        ("variants_after_ld_pruning", prune_in),
        ("pca_preprocessed_samples", final_samples),
        ("pca_preprocessed_variants", final_variants),
    ]

    out = Path(args.summary_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")


if __name__ == "__main__":
    main()
