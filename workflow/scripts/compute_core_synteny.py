#!/usr/bin/env python3
"""Compute per-contig core-synteny coverage from a self-vs-self minimap2
PAF alignment: the fraction of a contig's length covered by alignment
blocks to a *different* contig that is itself core-sized
(>= --core-min-length-bp).

High coverage means broad collinearity/duplication with a core chromosome
- evidence AGAINST the contig being a distinct accessory element (it is
more likely an assembly fragment/duplicate of a core chromosome). Low or
zero coverage is one of the criteria for a strong mini-chromosome candidate
(Starfish_und_MiniChromosomen_Analyseplan.md, C1/C2; user-supplied
criteria: "keine breite, kollineare Zuordnung zu Core-Chromosomen").

This does not detect redundancy between two non-core (e.g. both small)
contigs - that is a separate, not-yet-implemented fragmentation check.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def read_contig_ids(fasta: Path) -> list[str]:
    contig_ids = []
    with fasta.open("r") as infile:
        for line in infile:
            if line.startswith(">"):
                contig_ids.append(line[1:].split(maxsplit=1)[0].rstrip("\n"))
    return contig_ids


def _merge_and_sum(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals = sorted(intervals)
    covered = 0
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            covered += current_end - current_start
            current_start, current_end = start, end
    covered += current_end - current_start
    return covered


def compute_core_synteny(
    paf: Path, contig_ids: list[str], core_min_length_bp: int
) -> dict[str, float]:
    query_lengths: dict[str, int] = {}
    intervals_by_query: dict[str, list[tuple[int, int]]] = {cid: [] for cid in contig_ids}

    with paf.open("r") as infile:
        for line in infile:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            query_name, query_len, query_start, query_end = (
                fields[0],
                int(fields[1]),
                int(fields[2]),
                int(fields[3]),
            )
            target_name, target_len = fields[5], int(fields[6])

            query_lengths[query_name] = query_len

            if target_name == query_name:
                continue
            if target_len < core_min_length_bp:
                continue

            intervals_by_query.setdefault(query_name, []).append((query_start, query_end))

    coverage: dict[str, float] = {}
    for contig_id in contig_ids:
        query_len = query_lengths.get(contig_id)
        if query_len is None or query_len == 0:
            coverage[contig_id] = 0.0
            continue
        covered_bp = _merge_and_sum(intervals_by_query.get(contig_id, []))
        coverage[contig_id] = covered_bp / query_len

    return coverage


def write_coverage(coverage: dict[str, float], output: Path) -> None:
    with output.open("w") as outfile:
        outfile.write("contig_id\tcore_synteny_coverage\n")
        for contig_id, value in coverage.items():
            outfile.write(f"{contig_id}\t{value:.4f}\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paf", required=True, type=Path, help="Self-vs-self minimap2 PAF")
    parser.add_argument(
        "--fasta", required=True, type=Path, help="Normalized FASTA (for the full contig list)"
    )
    parser.add_argument("--core-min-length-bp", type=int, default=4_000_000)
    parser.add_argument("--output", required=True, type=Path, help="Output TSV path")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    contig_ids = read_contig_ids(args.fasta)
    coverage = compute_core_synteny(args.paf, contig_ids, args.core_min_length_bp)
    write_coverage(coverage, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
