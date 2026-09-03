#!/usr/bin/env python3
"""Section 7.3: classify every genomic window of every panel genome into
a region type, combining evidence already computed in earlier phases:

- gene-based prevalence (Orthogroups.tsv, Section 7.1)
- repeat density (repeats.smk, Section 6.3)
- contig-level mini-/accessory-chromosome signal (small + gene-poor,
  same logic used manually in docs/decisions.md 2026-09-02)
- starship_like candidate windows (starships.smk, Section 7.4)

Priority when several signals apply to the same window: starship_like
> accessory_chromosome (whole flagged contig) > subtelomeric_dynamic >
gene-prevalence-based core/shell/private > repeat_ambiguous >
unclassified. SyRI synteny (only available for 3/5 genomes, see
docs/decisions.md) is NOT yet folded in - documented as a follow-up.
"""
import argparse
import csv
from collections import defaultdict
from statistics import mean, median


def parse_fai(path):
    lengths = {}
    with open(path) as f:
        for line in f:
            cols = line.split("\t")
            lengths[cols[0]] = int(cols[1])
    return lengths


def parse_repeat_windows(path):
    windows = []
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            windows.append((cols[0], int(cols[1]), int(cols[2]), float(cols[-1])))
    return windows


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


def parse_orthogroups(path, genome_ids):
    """Return gene_id -> prevalence_count (how many of the 5 genomes carry
    its orthogroup), by parsing the Orthogroups.tsv gene-list columns."""
    gene_prevalence = {}
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        col_genomes = header[1:]
        for row in reader:
            cells = row[1:]
            present_genomes = sum(1 for c in cells if c.strip())
            for genome, cell in zip(col_genomes, cells):
                if not cell.strip():
                    continue
                for gid in cell.split(","):
                    gid = gid.strip()
                    if gid:
                        gene_prevalence[gid] = present_genomes
    return gene_prevalence


def mrna_to_gene_id(gff_path):
    """Liftoff/gffread output IDs genes as many mRNA-suffixed variants in
    OrthoFinder's proteome (e.g. 'rna-gnl|...|PoMZ_08913-RA_mrna'), while
    our GFF3 'gene' features use 'gene-PoMZ_08913'. Map gene-window
    overlap results (gene IDs) to the corresponding mRNA ID used in the
    proteome/orthogroups, via the locus_tag shared by both features."""
    locus_to_mrna = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "mRNA":
                continue
            attrs = dict(a.split("=", 1) for a in cols[8].split(";") if "=" in a)
            locus = attrs.get("locus_tag")
            mrna_id = attrs.get("ID")
            if locus and mrna_id:
                locus_to_mrna[locus] = mrna_id
    gene_to_locus = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attrs = dict(a.split("=", 1) for a in cols[8].split(";") if "=" in a)
            gid = attrs.get("ID")
            locus = attrs.get("locus_tag")
            if gid and locus:
                gene_to_locus[gid] = locus
    return {gid: locus_to_mrna.get(locus) for gid, locus in gene_to_locus.items()}


def parse_starship_windows(path):
    windows = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r["classification"] != "starship_like":
                continue
            windows[r["contig"]].append((int(r["window_start"]), int(r["window_end"])))
    return windows


def overlaps_any(intervals, start, end):
    return any(e > start and s < end for s, e in intervals)


def find_accessory_contigs(fai, genes_by_contig):
    lengths = list(fai.values())
    med_len = median(lengths)
    densities = {}
    for contig, length in fai.items():
        n_genes = len(genes_by_contig.get(contig, []))
        densities[contig] = n_genes / (length / 1_000_000) if length > 0 else 0
    mean_density = mean(densities.values()) if densities else 0
    accessory = set()
    for contig, length in fai.items():
        if length < 0.5 * med_len and densities[contig] < 0.5 * mean_density:
            accessory.add(contig)
    return accessory


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-ids", nargs="+", required=True)
    ap.add_argument("--fai-template", required=True)
    ap.add_argument("--repeat-windows-template", required=True)
    ap.add_argument("--repeat-per-contig-template", required=True)
    ap.add_argument("--gff-template", required=True)
    ap.add_argument("--orthogroups", required=True)
    ap.add_argument("--starship-candidates", required=True)
    ap.add_argument("--strict-core-fraction", type=float, required=True)
    ap.add_argument("--soft-core-fraction", type=float, required=True)
    ap.add_argument("--shell-min-fraction", type=float, required=True)
    ap.add_argument("--subtelomere-bp", type=int, default=50000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    n_genomes = len(args.genome_ids)
    gene_prevalence = parse_orthogroups(args.orthogroups, args.genome_ids)
    starship_windows = parse_starship_windows(args.starship_candidates)

    rows = []
    counts = defaultdict(int)

    for g in args.genome_ids:
        fai = parse_fai(args.fai_template.format(genome_id=g))
        rep_windows = parse_repeat_windows(args.repeat_windows_template.format(genome_id=g))
        genes_by_contig = parse_gff_genes(args.gff_template.format(genome_id=g))
        gene_to_mrna = mrna_to_gene_id(args.gff_template.format(genome_id=g))
        accessory_contigs = find_accessory_contigs(fai, genes_by_contig)

        with open(args.repeat_per_contig_template.format(genome_id=g)) as f:
            reader = csv.DictReader(f, delimiter="\t")
            fracs = [float(r["repeat_fraction"]) for r in reader]
        genome_avg_repeat = mean(fracs) if fracs else 0.0

        for contig, start, end, rep_frac in rep_windows:
            contig_len = fai.get(contig, end)
            genes_here = [
                (gs, ge, gid) for gs, ge, gid in genes_by_contig.get(contig, [])
                if ge > start and gs < end
            ]
            prevalences = []
            for gs, ge, gid in genes_here:
                mrna_id = gene_to_mrna.get(gid)
                if mrna_id and mrna_id in gene_prevalence:
                    prevalences.append(gene_prevalence[mrna_id])
            mean_prevalence_fraction = (mean(prevalences) / n_genomes) if prevalences else None

            near_telomere = start < args.subtelomere_bp or (contig_len - end) < args.subtelomere_bp
            repeat_enriched = rep_frac > genome_avg_repeat * 1.5

            if overlaps_any(starship_windows.get(contig, []), start, end):
                region_class = "starship_like"
            elif contig in accessory_contigs:
                region_class = "accessory_chromosome"
            elif near_telomere and repeat_enriched:
                region_class = "subtelomeric_dynamic"
            elif mean_prevalence_fraction is not None:
                if mean_prevalence_fraction >= args.strict_core_fraction:
                    region_class = "strict_core"
                elif mean_prevalence_fraction >= args.soft_core_fraction:
                    region_class = "soft_core"
                elif mean_prevalence_fraction >= args.shell_min_fraction:
                    region_class = "shell"
                else:
                    region_class = "private_accessory"
            elif repeat_enriched:
                region_class = "repeat_ambiguous"
            else:
                region_class = "unclassified"

            counts[region_class] += 1
            rows.append({
                "genome_id": g, "contig": contig, "start": start, "end": end,
                "repeat_fraction": f"{rep_frac:.4f}",
                "n_genes": len(genes_here),
                "mean_orthogroup_prevalence_fraction":
                    f"{mean_prevalence_fraction:.2f}" if mean_prevalence_fraction is not None else "NA",
                "region_class": region_class,
            })

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} windows to {args.out}")
    for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {k}: {v} ({100*v/len(rows):.1f}%)")


if __name__ == "__main__":
    main()
