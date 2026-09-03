#!/usr/bin/env python3
"""Section 13.1: assign each candidate panel region a stable candidate_id
and enrich with captain-gene/cargo-gene evidence where available.

Candidate classes (per document 13.1): accessory regions (ACC_), mini-
chromosomes (MCHR_), starship-like (STAR_), subtelomeric dynamic
(SUBTEL_). Our panel only carries three candidate classes after Section
7.3/8 (accessory_chromosome, starship_like, private_accessory) -
mini_chromosome/subtelomeric_dynamic were never assigned (no verified
telomere boundaries, see docs/decisions.md) so MCHR_/SUBTEL_ IDs do not
occur here; this is documented, not a bug.
"""
import argparse
import csv


def parse_source(source_interval):
    # e.g. "7015__CM113017.1:2030001-3040000"
    genome_contig, coords = source_interval.rsplit(":", 1)
    start, end = coords.split("-")
    return genome_contig, int(start), int(end)


def load_yr_hits(path):
    hits = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "mRNA":
                continue
            seqid, start, end = cols[0], int(cols[3]), int(cols[4])
            name = next((a.split("=", 1)[1] for a in cols[8].split(";") if a.startswith("Name=")), seqid)
            hits.append((seqid, start, end, name))
    return hits


def load_cargo_counts(path):
    counts = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            counts[r["yr_name"]] = int(r["n_cargo_genes"])
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--yr-gff", required=True)
    ap.add_argument("--starship-candidates", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    candidate_classes = {"accessory_chromosome", "starship_like", "private_accessory"}
    yr_hits = load_yr_hits(args.yr_gff)
    cargo_counts = load_cargo_counts(args.starship_candidates)

    rows = []
    with open(args.manifest) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["region_type"] not in candidate_classes:
                continue
            genome_contig, start, end = parse_source(r["source_interval"])

            captain_gene_id = "NA"
            cargo_gene_count = 0
            if r["region_type"] == "starship_like":
                for seqid, hstart, hend, name in yr_hits:
                    if seqid == genome_contig and hend > start and hstart < end:
                        captain_gene_id = name
                        cargo_gene_count = cargo_counts.get(name, 0)
                        break

            if r["region_type"] == "starship_like":
                annotation_confidence = "high" if captain_gene_id != "NA" and cargo_gene_count >= 1 else "medium"
            elif r["region_type"] in ("accessory_chromosome", "private_accessory"):
                annotation_confidence = "medium"
            else:
                annotation_confidence = "low"

            rows.append({
                "candidate_id": r["region_cluster"],
                "panel_id": r["panel_id"],
                "region_type": r["region_type"],
                "panel_start": 1,
                "panel_end": r["length_bp"],
                "length_bp": r["length_bp"],
                "representative_reference": r["representative_reference"],
                "reference_isolates": r["reference_isolates"],
                "n_reference_genomes": r["n_reference_genomes"],
                "host_groups": r["host_groups"],
                "core_accessory_class": r["region_type"],
                "repeat_fraction": r["repeat_fraction"],
                "captain_gene_id": captain_gene_id,
                "cargo_gene_count": cargo_gene_count,
                "annotation_confidence": annotation_confidence,
            })

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} candidate regions to {args.out}")
    from collections import Counter
    print(Counter(r["region_type"] for r in rows))


if __name__ == "__main__":
    main()
