#!/usr/bin/env python3
"""Calculate per-contig length, GC fraction and gene count as a first step
towards mini-chromosome (mChr) candidate metrics.

Further columns from the mChr metric table (core_gene_count,
repeat_fraction, TE_fraction, secreted_protein_count,
effector_candidate_count, median_depth, depth_ratio_to_core, mchr_score,
classification) require additional annotation/coverage inputs that are not
yet available for all isolates and are added once those tools are wired up.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

GC_BASES = set("GCgc")
ACGT_BASES = set("ACGTacgt")


def read_contig_lengths_and_gc(fasta: Path) -> dict[str, tuple[int, int, int]]:
    """Return {contig_id: (length_bp, gc_count, acgt_count)} in FASTA order."""
    metrics: dict[str, list[int]] = {}
    contig_id: str | None = None

    with fasta.open("r") as infile:
        for line in infile:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                contig_id = header.split(maxsplit=1)[0]
                metrics[contig_id] = [0, 0, 0]
                continue

            if contig_id is None:
                continue

            seq = line.strip()
            entry = metrics[contig_id]
            entry[0] += len(seq)
            entry[1] += sum(1 for base in seq if base in GC_BASES)
            entry[2] += sum(1 for base in seq if base in ACGT_BASES)

    return {contig_id: tuple(values) for contig_id, values in metrics.items()}


def read_gene_counts(gff: Path) -> dict[str, int]:
    """Return {contig_id: number of 'gene' features} from a GFF file."""
    counts: dict[str, int] = {}

    with gff.open("r") as infile:
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[2] != "gene":
                continue

            counts[fields[0]] = counts.get(fields[0], 0) + 1

    return counts


def read_assembly_units(sequence_report: Path, id_map: Path) -> dict[str, str]:
    """Return {new_id: assembly_unit} by joining the NCBI sequence report
    (keyed on the original accession, e.g. 'non-nuclear' for the
    mitochondrial genome) through the old_id -> new_id FASTA header mapping.
    """
    old_to_new: dict[str, str] = {}
    with id_map.open("r") as mapfile:
        next(mapfile)  # header
        for line in mapfile:
            if not line.strip():
                continue
            old_id, new_id = line.rstrip("\n").split("\t")
            old_to_new[old_id] = new_id

    units: dict[str, str] = {}
    with sequence_report.open("r") as infile:
        for line in infile:
            if not line.strip():
                continue
            record = json.loads(line)
            old_id = record.get("genbankAccession")
            new_id = old_to_new.get(old_id)
            if new_id is not None:
                units[new_id] = record.get("assemblyUnit", "")

    return units


def calculate_contig_metrics(
    fasta: Path,
    isolate_id: str,
    gff: Path | None = None,
    sequence_report: Path | None = None,
    id_map: Path | None = None,
) -> list[dict[str, object]]:
    contig_stats = read_contig_lengths_and_gc(fasta)
    gene_counts = read_gene_counts(gff) if gff is not None else None
    assembly_units = (
        read_assembly_units(sequence_report, id_map)
        if sequence_report is not None and id_map is not None
        else None
    )

    rows = []
    for contig_id, (length_bp, gc_count, acgt_count) in contig_stats.items():
        gc_fraction = gc_count / acgt_count if acgt_count > 0 else ""
        gene_count = gene_counts.get(contig_id, 0) if gene_counts is not None else ""
        assembly_unit = (
            assembly_units.get(contig_id, "") if assembly_units is not None else ""
        )
        rows.append(
            {
                "isolate_id": isolate_id,
                "contig_id": contig_id,
                "length_bp": length_bp,
                "gc_fraction": gc_fraction,
                "gene_count": gene_count,
                "assembly_unit": assembly_unit,
            }
        )
    return rows


def write_contig_metrics(rows: list[dict[str, object]], output: Path) -> None:
    columns = [
        "isolate_id",
        "contig_id",
        "length_bp",
        "gc_fraction",
        "gene_count",
        "assembly_unit",
    ]
    with output.open("w") as outfile:
        outfile.write("\t".join(columns) + "\n")
        for row in rows:
            values = [
                f"{row[col]:.4f}" if isinstance(row[col], float) else str(row[col])
                for col in columns
            ]
            outfile.write("\t".join(values) + "\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path, help="Normalized FASTA file")
    parser.add_argument(
        "--isolate-id", required=True, help="Isolate identifier, e.g. GCA_004346965.1"
    )
    parser.add_argument(
        "--gff",
        type=Path,
        default=None,
        help="Normalized GFF file (optional; enables gene_count)",
    )
    parser.add_argument(
        "--sequence-report",
        type=Path,
        default=None,
        help="Original NCBI sequence_report.jsonl (optional, used with --id-map "
        "to flag non-nuclear/mitochondrial contigs via assembly_unit)",
    )
    parser.add_argument(
        "--id-map",
        type=Path,
        default=None,
        help="old_id/new_id TSV from normalize_fasta_headers.py (required with "
        "--sequence-report)",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output TSV path")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    rows = calculate_contig_metrics(
        args.fasta, args.isolate_id, args.gff, args.sequence_report, args.id_map
    )
    write_contig_metrics(rows, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
