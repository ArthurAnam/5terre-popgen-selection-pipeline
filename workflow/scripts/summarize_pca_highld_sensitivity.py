#!/usr/bin/env python3

import argparse
from pathlib import Path


def count_nonempty(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


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
    parser.add_argument("--masked-snplist", required=True)
    parser.add_argument("--unmasked-prune-in", required=True)
    parser.add_argument("--masked-prune-in", required=True)
    parser.add_argument("--masked-prune-out", required=True)
    parser.add_argument("--final-pvar", required=True)
    parser.add_argument("--final-psam", required=True)
    parser.add_argument("--regions", required=True)
    parser.add_argument("--maf-threshold", type=float, required=True)
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
    masked_variants = count_nonempty(Path(args.masked_snplist))
    unmasked_prune_in = count_nonempty(Path(args.unmasked_prune_in))
    masked_prune_in = count_nonempty(Path(args.masked_prune_in))
    masked_prune_out = count_nonempty(Path(args.masked_prune_out))
    final_variants = count_pvar(Path(args.final_pvar))
    final_samples = count_psam(Path(args.final_psam))
    region_count = count_nonempty(Path(args.regions))

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
    if masked_variants > maf_variants:
        raise SystemExit("ERROR: high-LD masking increased the marker count")
    if masked_prune_in + masked_prune_out != masked_variants:
        raise SystemExit(
            "ERROR: masked prune.in + prune.out does not equal the masked "
            f"MAF-passing count ({masked_prune_in} + {masked_prune_out} != "
            f"{masked_variants})"
        )
    if final_variants != masked_prune_in:
        raise SystemExit(
            "ERROR: final masked PVAR count does not equal masked prune.in count "
            f"({final_variants} != {masked_prune_in})"
        )
    if final_samples != input_samples:
        raise SystemExit(
            "ERROR: high-LD sensitivity preprocessing unexpectedly changed the "
            f"sample count ({final_samples} != {input_samples})"
        )

    rows = [
        ("harmonized_input_variants", input_variants),
        ("harmonized_input_samples", input_samples),
        ("maf_threshold", args.maf_threshold),
        ("variants_after_maf", maf_variants),
        ("high_ld_region_count", region_count),
        ("variants_removed_by_high_ld_mask", maf_variants - masked_variants),
        ("variants_entering_masked_ld_pruning", masked_variants),
        ("ld_window_snps", args.window_snps),
        ("ld_step_snps", args.step_snps),
        ("ld_r2_threshold", args.r2_threshold),
        ("plink2_indep_order", args.indep_order),
        ("variants_removed_by_masked_ld_pruning", masked_prune_out),
        ("masked_variants_after_ld_pruning", masked_prune_in),
        ("unmasked_variants_after_ld_pruning", unmasked_prune_in),
        ("masked_minus_unmasked_final_variants", masked_prune_in - unmasked_prune_in),
        ("pca_preprocessed_samples_masked", final_samples),
        ("pca_preprocessed_variants_masked", final_variants),
    ]

    out = Path(args.summary_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")


if __name__ == "__main__":
    main()
