#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path

DNA = set("ACGT")

def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt", encoding="utf-8")

def read_pop(path, expected_pop):
    ids = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            q = line.split()
            if len(q) != 2:
                raise SystemExit(f"ERROR: malformed population row in {path}")
            if q[1] != expected_pop:
                raise SystemExit(f"ERROR: unexpected population label {q[1]} in {path}")
            ids.append(q[0])
    if not ids or len(ids) != len(set(ids)):
        raise SystemExit(f"ERROR: empty/duplicate population IDs in {path}")
    return ids

def expected_ids_from_bim(path):
    ids = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                q = line.split()
                if len(q) < 6:
                    raise SystemExit(f"ERROR: malformed BIM {path}")
                ids.append(q[1])
    return ids

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--population", required=True)
    p.add_argument("--expected-samples", type=int, required=True)
    p.add_argument("--pop-file", required=True)
    p.add_argument("--vcfs", nargs="+", required=True)
    p.add_argument("--maf-threshold", type=float, default=0.05)
    p.add_argument("--max-site-missing", type=float, default=0.05)
    p.add_argument("--require-maf", action="store_true")
    p.add_argument("--expected-bim")
    p.add_argument("--summary-out", required=True)
    p.add_argument("--by-chrom-out", required=True)
    p.add_argument("--panel-bim-out", required=True)
    args = p.parse_args()

    pop_ids = read_pop(args.pop_file, args.population)
    if len(pop_ids) != args.expected_samples:
        raise SystemExit(f"ERROR: {args.population}: expected {args.expected_samples} samples, found {len(pop_ids)}")

    rows = []
    panel_ids = []
    total_sites = total_missing_gt = total_called_gt = 0
    total_unphased = total_nonbinary = total_bad_alleles = 0
    total_duplicate = total_unsorted = total_maf_below = total_missing_above = 0
    sample_order = None

    Path(args.panel_bim_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.panel_bim_out, "w", encoding="utf-8") as bout:
        for path in args.vcfs:
            chrom = None
            n_sites = n_missing = n_called = n_unphased = n_nonbinary = 0
            n_bad_alleles = n_dup = n_unsorted = n_maf_below = n_missing_above = 0
            prev = None
            seen = set()

            with open_text(path) as f:
                header_samples = None
                for line in f:
                    if line.startswith("#CHROM"):
                        header_samples = line.rstrip("\n").split("\t")[9:]
                        break
                if header_samples is None:
                    raise SystemExit(f"ERROR: no #CHROM header in {path}")
                if header_samples != pop_ids:
                    raise SystemExit(f"ERROR: {args.population}: VCF sample order != population file in {path}")
                if sample_order is None:
                    sample_order = header_samples
                elif header_samples != sample_order:
                    raise SystemExit(f"ERROR: {args.population}: inconsistent sample order across chromosomes")

                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    q = line.rstrip("\n").split("\t")
                    if len(q) != 9 + args.expected_samples:
                        raise SystemExit(f"ERROR: malformed VCF row in {path}")
                    c, pos, ref, alt = q[0], int(q[1]), q[3], q[4]
                    if chrom is None:
                        chrom = c
                    elif c != chrom:
                        raise SystemExit(f"ERROR: multiple contigs in {path}")
                    if len(ref) != 1 or len(alt) != 1 or ref not in DNA or alt not in DNA or ref == alt:
                        n_bad_alleles += 1
                    key = (c, pos, ref, alt)
                    if key in seen:
                        n_dup += 1
                    seen.add(key)
                    if prev is not None and pos <= prev:
                        n_unsorted += 1
                    prev = pos

                    fmt = q[8].split(":")
                    if "GT" not in fmt:
                        raise SystemExit(f"ERROR: GT missing from FORMAT in {path}")
                    gi = fmt.index("GT")
                    ac = an = 0
                    missing_gt = 0
                    for sample in q[9:]:
                        vals = sample.split(":")
                        gt = vals[gi] if gi < len(vals) else "."
                        if "." in gt:
                            missing_gt += 1
                            n_missing += 1
                            continue
                        if "/" in gt:
                            n_unphased += 1
                        alleles = gt.replace("|", "/").split("/")
                        if len(alleles) != 2 or any(a not in ("0", "1") for a in alleles):
                            n_nonbinary += 1
                            continue
                        ac += sum(a == "1" for a in alleles)
                        an += 2
                        n_called += 1

                    miss_frac = missing_gt / args.expected_samples
                    if miss_frac > args.max_site_missing + 1e-12:
                        n_missing_above += 1
                    if an == 0:
                        n_maf_below += 1
                    else:
                        af = ac / an
                        maf = min(af, 1.0 - af)
                        if maf + 1e-12 < args.maf_threshold:
                            n_maf_below += 1

                    vid = f"{c}:{pos}:{ref}:{alt}"
                    panel_ids.append(vid)
                    bout.write(f"{c}\t{vid}\t0\t{pos}\t{alt}\t{ref}\n")
                    n_sites += 1

            if chrom is None:
                raise SystemExit(f"ERROR: no variants in {path}")
            rows.append({
                "chromosome": chrom,
                "sites": n_sites,
                "samples": args.expected_samples,
                "called_genotypes": n_called,
                "missing_genotypes": n_missing,
                "unphased_called_genotypes": n_unphased,
                "nonbinary_called_genotypes": n_nonbinary,
                "bad_snp_alleles": n_bad_alleles,
                "duplicate_sites": n_dup,
                "nonincreasing_positions": n_unsorted,
                "sites_maf_below_threshold": n_maf_below,
                "sites_missingness_above_threshold": n_missing_above,
            })
            total_sites += n_sites
            total_missing_gt += n_missing
            total_called_gt += n_called
            total_unphased += n_unphased
            total_nonbinary += n_nonbinary
            total_bad_alleles += n_bad_alleles
            total_duplicate += n_dup
            total_unsorted += n_unsorted
            total_maf_below += n_maf_below
            total_missing_above += n_missing_above

    expected_match = "NA"
    first_mismatch = "NA"
    if args.expected_bim:
        exp = expected_ids_from_bim(args.expected_bim)
        expected_match = "PASS" if exp == panel_ids else "FAIL"
        if exp != panel_ids:
            lim = min(len(exp), len(panel_ids))
            for i in range(lim):
                if exp[i] != panel_ids[i]:
                    first_mismatch = f"index={i};expected={exp[i]};observed={panel_ids[i]}"
                    break
            else:
                first_mismatch = f"length_expected={len(exp)};length_observed={len(panel_ids)}"

    hard_fail = any([
        total_unphased,
        total_nonbinary,
        total_bad_alleles,
        total_duplicate,
        total_unsorted,
        total_missing_above,
        args.require_maf and total_maf_below,
        expected_match == "FAIL",
    ])

    Path(args.by_chrom_out).parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with open(args.by_chrom_out, "w", encoding="utf-8", newline="") as out:
        w = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    with open(args.summary_out, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["audit_status", "FAIL" if hard_fail else "PASS"])
        w.writerow(["population", args.population])
        w.writerow(["autosomes", len(rows)])
        w.writerow(["samples", args.expected_samples])
        w.writerow(["total_sites", total_sites])
        w.writerow(["called_genotypes", total_called_gt])
        w.writerow(["missing_genotypes", total_missing_gt])
        w.writerow(["unphased_called_genotypes", total_unphased])
        w.writerow(["nonbinary_called_genotypes", total_nonbinary])
        w.writerow(["bad_snp_alleles", total_bad_alleles])
        w.writerow(["duplicate_sites", total_duplicate])
        w.writerow(["nonincreasing_positions", total_unsorted])
        w.writerow(["maf_threshold", args.maf_threshold])
        w.writerow(["sites_maf_below_threshold", total_maf_below])
        w.writerow(["max_site_missingness", args.max_site_missing])
        w.writerow(["sites_missingness_above_threshold", total_missing_above])
        w.writerow(["expected_ld_panel_exact_id_order", expected_match])
        w.writerow(["expected_ld_panel_first_mismatch", first_mismatch])

    if hard_fail:
        raise SystemExit(f"ERROR: {args.population} LASSI input audit failed; inspect {args.summary_out}")

if __name__ == "__main__":
    main()
