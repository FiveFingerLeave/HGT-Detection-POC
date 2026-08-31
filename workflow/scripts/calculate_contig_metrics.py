#!/usr/bin/env python3
"""Calculate per-contig length, GC fraction, gene count, repeat fraction and
telomere-completeness as a first step towards mini-chromosome (mChr)
candidate metrics.

`repeat_fraction` is the soft-masked base fraction from windowmasker
(genome-self-content based masking, no curated repeat library). This is a
coarse repeat-content proxy, not a curated TE-family annotation ("TE_fraction"
in the target schema) - that requires a tool like RepeatMasker/EDTA with a
species-appropriate repeat library and is deferred.

`telomere_start`/`telomere_end` flag whether a tandem run of the canonical
(TTAGGG)n fungal telomeric repeat (either orientation) is found within
`--telomere-window-bp` of each contig end - evidence that a long-read
assembly captured the true chromosome end there, per
Starfish_und_MiniChromosomen_Analyseplan.md.

Further columns from the mChr metric table (core_gene_count, TE_fraction,
secreted_protein_count, effector_candidate_count, median_depth,
depth_ratio_to_core, mchr_score, core-synteny) require additional
annotation/coverage/alignment inputs that are not yet available for all
isolates and are added once those tools are wired up.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

GC_BASES = set("GCgc")
ACGT_BASES = set("ACGTacgt")

# Canonical (TTAGGG)n telomeric repeat used by M. oryzae and most other
# filamentous ascomycetes, checked in both orientations since a contig's
# assembled strand/direction relative to the chromosome end is not known
# a priori.
TELOMERE_MOTIFS = ("TTAGGG", "CCCTAA")


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


def read_masked_fractions(masked_fasta: Path) -> dict[str, float]:
    """Return {contig_id: fraction of soft-masked (lowercase) bases} from a
    windowmasker-masked FASTA. Contig order/IDs must match the input FASTA.
    """
    fractions: dict[str, list[int]] = {}
    contig_id: str | None = None

    with masked_fasta.open("r") as infile:
        for line in infile:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                contig_id = header.split(maxsplit=1)[0]
                fractions[contig_id] = [0, 0]
                continue

            if contig_id is None:
                continue

            seq = line.strip()
            entry = fractions[contig_id]
            entry[0] += len(seq)
            entry[1] += sum(1 for base in seq if base.islower())

    return {
        contig_id: (masked / total if total > 0 else 0.0)
        for contig_id, (total, masked) in fractions.items()
    }


def _longest_tandem_run(seq: str, motif: str) -> int:
    """Longest run of immediately-consecutive (tandem) copies of motif in seq."""
    n = len(motif)
    best = current = 0
    i = 0
    while i <= len(seq) - n:
        if seq[i : i + n] == motif:
            current += 1
            best = max(best, current)
            i += n
        else:
            current = 0
            i += 1
    return best


def read_telomere_flags(
    fasta: Path, window_bp: int, min_repeats: int
) -> dict[str, tuple[bool, bool]]:
    """Return {contig_id: (telomere_start, telomere_end)}: whether a tandem
    run of >= min_repeats telomeric repeat units (in either orientation) is
    found within the first/last window_bp bases of the contig.
    """
    start_buffers: dict[str, str] = {}
    end_buffers: dict[str, str] = {}
    contig_id: str | None = None

    with fasta.open("r") as infile:
        for line in infile:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                contig_id = header.split(maxsplit=1)[0]
                start_buffers[contig_id] = ""
                end_buffers[contig_id] = ""
                continue

            if contig_id is None:
                continue

            seq = line.strip().upper()
            if len(start_buffers[contig_id]) < window_bp:
                start_buffers[contig_id] = (start_buffers[contig_id] + seq)[:window_bp]
            end_buffers[contig_id] = (end_buffers[contig_id] + seq)[-window_bp:]

    flags: dict[str, tuple[bool, bool]] = {}
    for cid in start_buffers:
        start_hit = any(
            _longest_tandem_run(start_buffers[cid], motif) >= min_repeats
            for motif in TELOMERE_MOTIFS
        )
        end_hit = any(
            _longest_tandem_run(end_buffers[cid], motif) >= min_repeats
            for motif in TELOMERE_MOTIFS
        )
        flags[cid] = (start_hit, end_hit)

    return flags


def read_sequence_report_fields(
    sequence_report: Path, id_map: Path
) -> dict[str, dict[str, str]]:
    """Return {new_id: {"assembly_unit": ..., "role": ...}} by joining the
    NCBI sequence report (keyed on the original accession) through the
    old_id -> new_id FASTA header mapping.

    `assembly_unit` flags non-nuclear (e.g. mitochondrial) sequences.
    `role` distinguishes 'assembled-molecule' (a full chromosome) from
    'unplaced-scaffold' (sequence the assembler could not anchor to a
    chromosome, often because it is highly repetitive) - a small, repeat-rich
    unplaced scaffold can look like an mChr candidate by the numbers but may
    simply be unassembled/unanchored repetitive sequence.
    """
    old_to_new: dict[str, str] = {}
    with id_map.open("r") as mapfile:
        next(mapfile)  # header
        for line in mapfile:
            if not line.strip():
                continue
            old_id, new_id = line.rstrip("\n").split("\t")
            old_to_new[old_id] = new_id

    fields: dict[str, dict[str, str]] = {}
    with sequence_report.open("r") as infile:
        for line in infile:
            if not line.strip():
                continue
            record = json.loads(line)
            old_id = record.get("genbankAccession")
            new_id = old_to_new.get(old_id)
            if new_id is not None:
                fields[new_id] = {
                    "assembly_unit": record.get("assemblyUnit", ""),
                    "role": record.get("role", ""),
                }

    return fields


def calculate_contig_metrics(
    fasta: Path,
    isolate_id: str,
    gff: Path | None = None,
    sequence_report: Path | None = None,
    id_map: Path | None = None,
    masked_fasta: Path | None = None,
    telomere_window_bp: int = 1000,
    telomere_min_repeats: int = 5,
) -> list[dict[str, object]]:
    contig_stats = read_contig_lengths_and_gc(fasta)
    gene_counts = read_gene_counts(gff) if gff is not None else None
    sequence_report_fields = (
        read_sequence_report_fields(sequence_report, id_map)
        if sequence_report is not None and id_map is not None
        else None
    )
    masked_fractions = (
        read_masked_fractions(masked_fasta) if masked_fasta is not None else None
    )
    telomere_flags = read_telomere_flags(fasta, telomere_window_bp, telomere_min_repeats)

    rows = []
    for contig_id, (length_bp, gc_count, acgt_count) in contig_stats.items():
        gc_fraction = gc_count / acgt_count if acgt_count > 0 else ""
        gene_count = gene_counts.get(contig_id, 0) if gene_counts is not None else ""
        report_fields = (
            sequence_report_fields.get(contig_id, {})
            if sequence_report_fields is not None
            else {}
        )
        repeat_fraction = (
            masked_fractions.get(contig_id, "") if masked_fractions is not None else ""
        )
        telomere_start, telomere_end = telomere_flags.get(contig_id, (False, False))
        rows.append(
            {
                "isolate_id": isolate_id,
                "contig_id": contig_id,
                "length_bp": length_bp,
                "gc_fraction": gc_fraction,
                "gene_count": gene_count,
                "repeat_fraction": repeat_fraction,
                "assembly_unit": report_fields.get("assembly_unit", ""),
                "role": report_fields.get("role", ""),
                "telomere_start": telomere_start,
                "telomere_end": telomere_end,
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
        "repeat_fraction",
        "assembly_unit",
        "role",
        "telomere_start",
        "telomere_end",
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
    parser.add_argument(
        "--masked-fasta",
        type=Path,
        default=None,
        help="windowmasker soft-masked version of --fasta (optional; enables "
        "repeat_fraction)",
    )
    parser.add_argument(
        "--telomere-window-bp",
        type=int,
        default=1000,
        help="Bases from each contig end to scan for telomeric repeats (default: 1000)",
    )
    parser.add_argument(
        "--telomere-min-repeats",
        type=int,
        default=5,
        help="Minimum tandem (TTAGGG)n/(CCCTAA)n repeat units to call a telomere "
        "(default: 5; real M. oryzae telomeres in test data showed 18-31 vs. "
        "0-1 for non-telomeric ends)",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output TSV path")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    rows = calculate_contig_metrics(
        args.fasta,
        args.isolate_id,
        args.gff,
        args.sequence_report,
        args.id_map,
        args.masked_fasta,
        args.telomere_window_bp,
        args.telomere_min_repeats,
    )
    write_contig_metrics(rows, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
