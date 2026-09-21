#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path


def read_psam_ids(path: Path):
    ids = []
    with path.open("r", encoding="utf-8") as handle:
        header = None
        iid_index = None
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            fields = line.split()
            if header is None:
                header = fields
                normalized = [x.lstrip("#") for x in header]
                if "IID" not in normalized:
                    raise SystemExit("ERROR: PSAM header does not contain IID")
                iid_index = normalized.index("IID")
                continue
            ids.append(fields[iid_index])
    if len(ids) != len(set(ids)):
        raise SystemExit("ERROR: duplicate sample IDs in PSAM")
    return ids


def read_eur_metadata(path: Path):
    mapping = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ID", "POP"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit("ERROR: EUR metadata must contain ID and POP columns")
        for row in reader:
            sample = row["ID"]
            pop = row["POP"]
            if sample in mapping:
                raise SystemExit(f"ERROR: duplicate EUR metadata ID: {sample}")
            mapping[sample] = pop
    return mapping


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psam", required=True)
    parser.add_argument("--eur-metadata", required=True)
    parser.add_argument("--expected-total", type=int, required=True)
    parser.add_argument("--expected-eur", type=int, required=True)
    parser.add_argument("--expected-ct", type=int, required=True)
    parser.add_argument("--annotations-out", required=True)
    parser.add_argument("--counts-out", required=True)
    args = parser.parse_args()

    sample_ids = read_psam_ids(Path(args.psam))
    eur = read_eur_metadata(Path(args.eur_metadata))

    if len(sample_ids) != args.expected_total:
        raise SystemExit(
            f"ERROR: expected {args.expected_total} joint PCA samples, found {len(sample_ids)}"
        )

    present_eur = [sid for sid in sample_ids if sid in eur]
    missing_eur = sorted(set(eur) - set(sample_ids))
    if missing_eur:
        raise SystemExit(
            f"ERROR: {len(missing_eur)} expected EUR metadata samples are absent from PSAM"
        )
    if len(present_eur) != args.expected_eur:
        raise SystemExit(
            f"ERROR: expected {args.expected_eur} EUR samples, found {len(present_eur)}"
        )

    rows = []
    for sid in sample_ids:
        if sid in eur:
            rows.append((sid, eur[sid], "1000G_EUR"))
        else:
            rows.append((sid, "CT", "Cinque_Terre"))

    counts = Counter(pop for _, pop, _ in rows)
    if counts["CT"] != args.expected_ct:
        raise SystemExit(
            f"ERROR: expected {args.expected_ct} Cinque Terre samples, found {counts['CT']}"
        )

    allowed = {"CT", "CEU", "FIN", "GBR", "IBS", "TSI"}
    unexpected = sorted(set(counts) - allowed)
    if unexpected:
        raise SystemExit(f"ERROR: unexpected PCA population labels: {unexpected}")

    annotations_out = Path(args.annotations_out)
    counts_out = Path(args.counts_out)
    annotations_out.parent.mkdir(parents=True, exist_ok=True)

    with annotations_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "population", "source"])
        writer.writerows(rows)

    order = ["CT", "CEU", "FIN", "GBR", "IBS", "TSI"]
    with counts_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["population", "n_samples"])
        for pop in order:
            writer.writerow([pop, counts[pop]])
        writer.writerow(["TOTAL", len(rows)])


if __name__ == "__main__":
    main()
