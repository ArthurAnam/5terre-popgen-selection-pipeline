#!/usr/bin/env python3
import argparse
import gzip
import re
from pathlib import Path

DNA = set("ACGT")
COMP = dict(zip("ATCG", "TAGC"))

def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt", encoding="utf-8")

def canon(pos, a, b):
    if a not in DNA or b not in DNA or a == b:
        raise SystemExit(f"ERROR: invalid SNP alleles at {pos}: {a}/{b}")
    direct = tuple(sorted((a, b)))
    comp = tuple(sorted((COMP[a], COMP[b])))
    return (int(pos), *min(direct, comp))

def sample_ids(path):
    rows = [x.split() for x in Path(path).read_text().splitlines() if x.strip()]
    if len(rows) < 3 or rows[0][:2] != ["ID_1", "ID_2"] or rows[1][:2] != ["0", "0"]:
        raise SystemExit(f"ERROR: malformed SHAPEIT2 sample file: {path}")
    if any(r[0] != r[1] for r in rows[2:]):
        raise SystemExit(f"ERROR: ID_1/ID_2 mismatch in {path}")
    ids = [r[0] for r in rows[2:]]
    if len(ids) != len(set(ids)):
        raise SystemExit("ERROR: duplicate sample IDs")
    return ids

def source_header_and_iterator(path):
    handle = open_text(path)
    samples = None
    def records():
        nonlocal samples
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
        if samples is None:
            raise SystemExit(f"ERROR: no #CHROM header in {path}")
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise SystemExit("ERROR: malformed source VCF row")
            yield fields
    it = records()
    # Prime the generator so the header is parsed before returning sample IDs.
    try:
        first = next(it)
    except StopIteration:
        raise SystemExit(f"ERROR: no variants in {path}")
    def chained():
        yield first
        yield from it
    return handle, lambda: samples, chained()

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--haps", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--source-vcf", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    ids = sample_ids(args.sample)
    source_handle, get_source_samples, source_iter = source_header_and_iterator(args.source_vcf)
    source_samples = get_source_samples()
    if source_samples != ids:
        raise SystemExit("ERROR: SHAPEIT2 sample order does not match source VCF")

    src = next(source_iter, None)
    n = 0
    prev_pos = None

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.haps) as hin, gzip.open(args.output, "wt", compresslevel=6) as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=SHAPEIT2_v2.r904_to_lassip_converter\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")

        for line in hin:
            if not line.strip():
                continue
            h = line.split()
            if len(h) != 5 + 2 * len(ids):
                raise SystemExit(f"ERROR: malformed HAPS row {n+1}")
            hchrom, hpos, ha, hb = h[0], int(h[2]), h[3], h[4]
            hk = canon(hpos, ha, hb)

            while src is not None and canon(src[1], src[3], src[4]) != hk:
                if int(src[1]) > hpos:
                    raise SystemExit(f"ERROR: phased site {hchrom}:{hpos} not found in source VCF")
                src = next(source_iter, None)

            if src is None:
                raise SystemExit(f"ERROR: source VCF ended before phased site {hchrom}:{hpos}")

            schrom, spos, ref, alt = src[0], int(src[1]), src[3], src[4]
            if spos != hpos:
                raise SystemExit(f"ERROR: position mismatch at phased site {hchrom}:{hpos}")

            same = (ha, hb) == (ref, alt) or (COMP[ha], COMP[hb]) == (ref, alt)
            rev = (hb, ha) == (ref, alt) or (COMP[hb], COMP[ha]) == (ref, alt)
            if same == rev:
                raise SystemExit(f"ERROR: ambiguous/unresolved allele orientation at {schrom}:{spos}")

            bits = h[5:]
            if any(x not in ("0", "1") for x in bits):
                raise SystemExit(f"ERROR: non-binary HAPS allele at {schrom}:{spos}")
            if rev:
                bits = ["1" if x == "0" else "0" for x in bits]

            if prev_pos is not None and hpos <= prev_pos:
                raise SystemExit(f"ERROR: non-increasing HAPS position at {schrom}:{spos}")
            prev_pos = hpos

            vid = f"{schrom}:{spos}:{ref}:{alt}"
            gts = [bits[2*i] + "|" + bits[2*i+1] for i in range(len(ids))]
            out.write(f"{schrom}\t{spos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
            n += 1
            src = next(source_iter, None)

    source_handle.close()
    if n == 0:
        raise SystemExit("ERROR: no phased variants written")

if __name__ == "__main__":
    main()
