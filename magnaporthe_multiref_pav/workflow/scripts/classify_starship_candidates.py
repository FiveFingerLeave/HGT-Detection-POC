#!/usr/bin/env python3
"""Section 7.4: synthesize YR/Captain hits (starfish annotate) with
window size, repeat context (repeats.smk), and cargo genes (Liftoff
GFF3, annotation.smk) into a conservative starship_like classification.

Deliberately does NOT require starfish's own internal name-lifting
(which failed silently on our Liftoff GFF3 - see docs/decisions.md) -
cargo genes are looked up directly against data/annotations/*.gff3 by
genomic overlap instead.
"""
import argparse
import csv
from collections import defaultdict


def parse_yr_gff(path):
    hits = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "mRNA":
                continue
            seqid, start, end, strand, attrs = cols[0], int(cols[3]), int(cols[4]), cols[6], cols[8]
            name = next((a.split("=", 1)[1] for a in attrs.split(";") if a.startswith("Name=")), seqid)
            genome_id = seqid.split("__", 1)[0]
            hits.append({"genome_id": genome_id, "contig": seqid, "start": start, "end": end,
                         "strand": strand, "yr_name": name})
    return hits


def parse_fai(path):
    lengths = {}
    with open(path) as f:
        for line in f:
            cols = line.split("\t")
            lengths[cols[0]] = int(cols[1])
    return lengths


def parse_repeat_windows(path):
    # bedtools coverage output: chrom start end n_overlaps bases_covered length fraction
    windows = defaultdict(list)
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            chrom, start, end, frac = cols[0], int(cols[1]), int(cols[2]), float(cols[-1])
            windows[chrom].append((start, end, frac))
    return windows


def mean_repeat_fraction(windows_for_contig, win_start, win_end):
    overlapping = [(s, e, f) for s, e, f in windows_for_contig if e > win_start and s < win_end]
    if not overlapping:
        return None
    total_overlap = sum(min(e, win_end) - max(s, win_start) for s, e, f in overlapping)
    if total_overlap <= 0:
        return None
    weighted = sum((min(e, win_end) - max(s, win_start)) * f for s, e, f in overlapping)
    return weighted / total_overlap


def parse_gff_genes(path):
    genes = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            seqid, start, end = cols[0], int(cols[3]), int(cols[4])
            gene_id = next((a.split("=", 1)[1] for a in cols[8].split(";") if a.startswith("ID=")), None)
            genes[seqid].append((start, end, gene_id))
    return genes


def genes_in_window(genes_for_contig, win_start, win_end):
    return [gid for s, e, gid in genes_for_contig if e > win_start and s < win_end]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--yr-gff", required=True)
    ap.add_argument("--genome-ids", nargs="+", required=True)
    ap.add_argument("--fai-template", required=True, help="e.g. data/references/{genome_id}.fa.fai")
    ap.add_argument("--repeat-windows-template", required=True,
                     help="e.g. results/repeats/{genome_id}_repeat_windows.bed")
    ap.add_argument("--gff-template", required=True, help="e.g. data/annotations/{genome_id}.gff3")
    ap.add_argument("--repeat-per-contig-template", required=True,
                     help="e.g. results/repeats/{genome_id}_repeat_per_contig.tsv")
    ap.add_argument("--min-region-bp", type=int, default=20000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    hits = parse_yr_gff(args.yr_gff)

    genome_avg_repeat = {}
    for g in args.genome_ids:
        with open(args.repeat_per_contig_template.format(genome_id=g)) as f:
            reader = csv.DictReader(f, delimiter="\t")
            fracs = [float(r["repeat_fraction"]) for r in reader]
        genome_avg_repeat[g] = sum(fracs) / len(fracs) if fracs else 0.0

    rows = []
    for g in args.genome_ids:
        fai = parse_fai(args.fai_template.format(genome_id=g))
        rep_windows = parse_repeat_windows(args.repeat_windows_template.format(genome_id=g))
        genes = parse_gff_genes(args.gff_template.format(genome_id=g))
        genome_hits = [h for h in hits if h["genome_id"] == g]
        for h in genome_hits:
            contig_len = fai.get(h["contig"], h["end"])
            hit_len = h["end"] - h["start"]
            pad = max(0, (args.min_region_bp - hit_len) // 2)
            win_start = max(0, h["start"] - pad)
            win_end = min(contig_len, h["end"] + pad)
            region_len = win_end - win_start

            rep_frac = mean_repeat_fraction(rep_windows.get(h["contig"], []), win_start, win_end)
            cargo_genes = genes_in_window(genes.get(h["contig"], []), win_start, win_end)
            # exclude the YR gene's own overlapping gene model, if any, from cargo count
            n_cargo = len(cargo_genes)

            repeat_enriched = rep_frac is not None and rep_frac > genome_avg_repeat[g] * 1.5
            meets_size = region_len >= args.min_region_bp
            has_cargo = n_cargo >= 1

            if meets_size and has_cargo and repeat_enriched:
                classification = "starship_like"
            elif has_cargo or repeat_enriched:
                classification = "duf3435_candidate_contextual"
            else:
                classification = "duf3435_candidate_only"

            rows.append({
                "genome_id": g,
                "contig": h["contig"],
                "yr_name": h["yr_name"],
                "yr_start": h["start"],
                "yr_end": h["end"],
                "strand": h["strand"],
                "window_start": win_start,
                "window_end": win_end,
                "window_bp": region_len,
                "meets_min_size_20kb": meets_size,
                "repeat_fraction_window": f"{rep_frac:.4f}" if rep_frac is not None else "NA",
                "genome_avg_repeat_fraction": f"{genome_avg_repeat[g]:.4f}",
                "repeat_enriched_1.5x": repeat_enriched,
                "n_cargo_genes": n_cargo,
                "classification": classification,
            })

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} candidates to {args.out}")
    counts = defaultdict(int)
    for r in rows:
        counts[r["classification"]] += 1
    for k, v in sorted(counts.items()):
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
