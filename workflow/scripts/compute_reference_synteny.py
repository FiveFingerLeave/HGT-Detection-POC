#!/usr/bin/env python3
"""Compute reference-panel synteny consensus for each contig of an isolate.

For each structural-tier reference genome, a minimap2 PAF (isolate as
query, reference as target) is intersected with that reference's core
contig list (see list_reference_core_contigs.py) to get a per-contig
coverage fraction against that reference. A contig "hits" a reference if
that coverage is at least --min-coverage-per-reference. `reference_core_hits`
is the number of references hit; `reference_core_consensus` is whether
that meets --consensus-min-hits.

Coverage against a reference is the BEST coverage to any SINGLE target
contig, not the sum across all of a reference's core contigs, and only
alignment blocks with mapping quality >= --min-mapq count. Without this, a
repeat/TE family shared with the core genome but scattered across many
different chromosomes at low, ambiguous mapping quality can accumulate
enough total query coverage to look like broad, confident collinearity
with "the core genome" even though no single alignment block represents
genuine 1:1 synteny to any one chromosome. (Observed on B71's own
mini-chromosome contig during development - see docs/decisions.md.)

This is a *different* question from compute_core_synteny.py's
core_synteny_coverage (which checks duplication within the SAME isolate's
own assembly): here, matching many independent reference genomes is
evidence the contig belongs to the conserved, species-wide core genome,
regardless of how it happens to look internally in this one isolate.
Failing to match any reference is one of the strongest available signals
for lineage-specific/accessory content.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def read_contig_ids_and_lengths(fasta: Path) -> dict[str, int]:
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


def compute_reference_coverage(
    paf: Path, core_contigs: set[str], min_mapq: int
) -> dict[str, dict[str, list[tuple[int, int]]]]:
    """Return {query_contig_id: {target_contig_id: [(start, end), ...]}},
    restricted to alignments whose target is one of this reference's core
    contigs and whose mapping quality is >= min_mapq."""
    intervals: dict[str, dict[str, list[tuple[int, int]]]] = {}
    with paf.open("r") as infile:
        for line in infile:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query_name, query_start, query_end = fields[0], int(fields[2]), int(fields[3])
            target_name = fields[5]
            mapq = int(fields[11])
            if target_name not in core_contigs or mapq < min_mapq:
                continue
            intervals.setdefault(query_name, {}).setdefault(target_name, []).append(
                (query_start, query_end)
            )
    return intervals


def compute_reference_synteny(
    query_fasta: Path,
    reference_pafs: list[Path],
    reference_core_contig_files: list[Path],
    min_coverage_per_reference: float,
    consensus_min_hits: int,
    min_mapq: int = 30,
) -> dict[str, dict[str, object]]:
    query_lengths = read_contig_ids_and_lengths(query_fasta)

    hits_per_contig: dict[str, int] = {cid: 0 for cid in query_lengths}

    for paf, core_contigs_file in zip(reference_pafs, reference_core_contig_files):
        core_contigs = set(core_contigs_file.read_text().split())
        intervals = compute_reference_coverage(paf, core_contigs, min_mapq)
        for contig_id, query_len in query_lengths.items():
            if query_len == 0:
                continue
            per_target = intervals.get(contig_id, {})
            best_coverage = max(
                (_merge_and_sum(ivs) / query_len for ivs in per_target.values()),
                default=0.0,
            )
            if best_coverage >= min_coverage_per_reference:
                hits_per_contig[contig_id] += 1

    return {
        contig_id: {
            "reference_core_hits": hits,
            "reference_core_consensus": hits >= consensus_min_hits,
        }
        for contig_id, hits in hits_per_contig.items()
    }


def write_reference_synteny(result: dict[str, dict[str, object]], output: Path) -> None:
    with output.open("w") as outfile:
        outfile.write("contig_id\treference_core_hits\treference_core_consensus\n")
        for contig_id, values in result.items():
            outfile.write(
                f"{contig_id}\t{values['reference_core_hits']}\t{values['reference_core_consensus']}\n"
            )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-fasta", required=True, type=Path)
    parser.add_argument(
        "--reference-paf", action="append", required=True, type=Path, dest="reference_pafs"
    )
    parser.add_argument(
        "--reference-core-contigs",
        action="append",
        required=True,
        type=Path,
        dest="reference_core_contig_files",
    )
    parser.add_argument("--min-coverage-per-reference", type=float, default=0.3)
    parser.add_argument("--consensus-min-hits", type=int, default=2)
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=30,
        help="Minimum minimap2 mapping quality for an alignment block to count "
        "(default: 30; excludes ambiguous/repetitive multi-mapping alignments)",
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    if len(args.reference_pafs) != len(args.reference_core_contig_files):
        raise SystemExit(
            "--reference-paf and --reference-core-contigs must be given the same "
            "number of times, in matching order"
        )

    result = compute_reference_synteny(
        args.query_fasta,
        args.reference_pafs,
        args.reference_core_contig_files,
        args.min_coverage_per_reference,
        args.consensus_min_hits,
        args.min_mapq,
    )
    write_reference_synteny(result, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
