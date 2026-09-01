#!/usr/bin/env python3
"""Classify each candidate region (Starship/core, from candidate_regions.bed)
by whether it is syntenic across a pair of reference genomes, using a
whole-genome minimap2 alignment (PAF) of the region's own genome (query)
against another reference genome (target).

Rationale: a horizontally-acquired mobile element (Starship) is expected to
show the classic insertion-polymorphism signature - its flanking sequence
aligns cleanly to the other genome (the insertion site is conserved), but
the region itself does not (the cargo is absent there). A region that
aligns cleanly end-to-end in the other genome too is not lineage-specific
in this pairwise comparison (shared ancestrally, or independently present in
both). A region where neither the region nor its flanks align is
uninformative (near a contig end, in a repeat-rich/divergent area) rather
than a synteny breakpoint.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def read_contig_names(fasta_path: Path) -> set[str]:
    names = set()
    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                names.add(line[1:].split()[0])
    return names


def read_regions(bed_path: Path, valid_contigs: set[str] | None = None) -> list[dict]:
    regions = []
    with open(bed_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            chrom, start, end, region_id, region_class = line.split("\t")
            if valid_contigs is not None and chrom not in valid_contigs:
                continue
            regions.append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "region_id": region_id,
                    "region_class": region_class,
                }
            )
    return regions


def parse_paf(paf_path: Path) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = {}
    with open(paf_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            qname, qstart, qend = fields[0], int(fields[2]), int(fields[3])
            intervals.setdefault(qname, []).append((qstart, qend))
    return intervals


def merge_intervals(intervals: list[tuple[int, int]]) -> list[list[int]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def aligned_fraction(intervals: list[list[int]], start: int, end: int) -> float:
    span = end - start
    if span <= 0:
        return 0.0
    covered = 0
    for istart, iend in intervals:
        overlap = min(iend, end) - max(istart, start)
        if overlap > 0:
            covered += overlap
    return covered / span


def classify_synteny(
    region_frac: float,
    upstream_frac: float,
    downstream_frac: float,
    present_threshold: float = 0.7,
    absent_threshold: float = 0.1,
    flank_threshold: float = 0.5,
) -> str:
    """See module docstring for the biological rationale of each label."""
    if region_frac >= present_threshold:
        return "syntenic_present"
    if (
        region_frac <= absent_threshold
        and upstream_frac >= flank_threshold
        and downstream_frac >= flank_threshold
    ):
        return "syntenic_empty_site"
    if upstream_frac < 0.3 and downstream_frac < 0.3:
        return "non_syntenic_region"
    return "ambiguous"


def evaluate_regions(
    regions: list[dict],
    paf_intervals: dict[str, list[tuple[int, int]]],
    flank_bp: int,
    genome: str,
    other_genome: str,
) -> list[dict]:
    rows = []
    for region in regions:
        chrom_intervals = merge_intervals(paf_intervals.get(region["chrom"], []))
        region_frac = aligned_fraction(chrom_intervals, region["start"], region["end"])
        up_frac = aligned_fraction(
            chrom_intervals, max(0, region["start"] - flank_bp), region["start"]
        )
        down_frac = aligned_fraction(chrom_intervals, region["end"], region["end"] + flank_bp)
        classification = classify_synteny(region_frac, up_frac, down_frac)
        rows.append(
            {
                "region_id": region["region_id"],
                "region_class": region["region_class"],
                "genome": genome,
                "other_genome": other_genome,
                "region_length_bp": region["end"] - region["start"],
                "region_aligned_fraction": round(region_frac, 4),
                "upstream_flank_aligned_fraction": round(up_frac, 4),
                "downstream_flank_aligned_fraction": round(down_frac, 4),
                "classification": classification,
            }
        )
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", required=True, type=Path, help="candidate_regions.bed")
    parser.add_argument(
        "--fasta", required=True, type=Path, help="Query genome FASTA (defines valid contigs)"
    )
    parser.add_argument("--paf", required=True, type=Path, help="minimap2 PAF: query vs other genome")
    parser.add_argument("--genome", required=True, help="Query genome/reference ID")
    parser.add_argument("--other-genome", required=True, help="Target genome/reference ID")
    parser.add_argument("--flank-bp", type=int, default=2000)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    valid_contigs = read_contig_names(args.fasta)
    regions = read_regions(args.regions, valid_contigs)
    paf_intervals = parse_paf(args.paf)
    rows = evaluate_regions(regions, paf_intervals, args.flank_bp, args.genome, args.other_genome)

    fieldnames = [
        "region_id",
        "region_class",
        "genome",
        "other_genome",
        "region_length_bp",
        "region_aligned_fraction",
        "upstream_flank_aligned_fraction",
        "downstream_flank_aligned_fraction",
        "classification",
    ]
    import csv

    with args.output.open("w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
