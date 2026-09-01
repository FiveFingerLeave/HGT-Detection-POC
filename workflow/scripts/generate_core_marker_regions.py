#!/usr/bin/env python3
"""Generate core-genome marker windows for one reference genome: evenly
spaced, fixed-size windows on its single largest contig (virtually
certain to be a core chromosome, not an accessory element), skipping any
window that overlaps a known candidate (Starship/mChr) region.

These are a coarse, tool-free stand-in for real single-copy conserved
markers (e.g. BUSCO genes) - see docs/decisions.md for why (BUSCO could
not be run in this environment: OOM). They serve the same purpose in
Section 3's core-vs-candidate clustering comparison: a background of loci
expected to follow vertical inheritance/host lineage cleanly.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def read_contig_lengths(fasta: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    contig_id: str | None = None
    with fasta.open("r") as infile:
        for line in infile:
            if line.startswith(">"):
                contig_id = line[1:].split(maxsplit=1)[0].rstrip("\n")
                lengths[contig_id] = 0
                continue
            if contig_id is not None:
                lengths[contig_id] += len(line.strip())
    return lengths


def read_excluded_regions(bed: Path) -> dict[str, list[tuple[int, int]]]:
    excluded: dict[str, list[tuple[int, int]]] = {}
    if not bed.exists():
        return excluded
    with bed.open("r") as infile:
        for line in infile:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            excluded.setdefault(chrom, []).append((start, end))
    return excluded


def _overlaps(start: int, end: int, intervals: list[tuple[int, int]]) -> bool:
    return any(start < iv_end and end > iv_start for iv_start, iv_end in intervals)


def generate_core_windows(
    contig_lengths: dict[str, int],
    excluded_regions: dict[str, list[tuple[int, int]]],
    genome_id: str,
    window_size: int,
    n_windows: int,
) -> list[tuple[str, int, int, str, str]]:
    if not contig_lengths:
        return []

    largest_contig = max(contig_lengths, key=contig_lengths.get)
    length = contig_lengths[largest_contig]
    excluded = excluded_regions.get(largest_contig, [])

    stride = max(window_size, length // (n_windows + 1))
    regions = []
    position = stride
    while len(regions) < n_windows and position + window_size <= length:
        start, end = position, position + window_size
        if not _overlaps(start, end, excluded):
            region_id = f"core_{genome_id}_{len(regions) + 1}"
            regions.append((largest_contig, start, end, region_id, "core"))
        position += stride

    return regions


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--genome-id", required=True)
    parser.add_argument("--exclude-bed", type=Path, default=None)
    parser.add_argument("--window-size", type=int, default=50_000)
    parser.add_argument("--n-windows", type=int, default=20)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    contig_lengths = read_contig_lengths(args.fasta)
    excluded = read_excluded_regions(args.exclude_bed) if args.exclude_bed else {}
    regions = generate_core_windows(
        contig_lengths, excluded, args.genome_id, args.window_size, args.n_windows
    )

    with args.output.open("w") as outfile:
        for chrom, start, end, region_id, region_class in regions:
            outfile.write(f"{chrom}\t{start}\t{end}\t{region_id}\t{region_class}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
