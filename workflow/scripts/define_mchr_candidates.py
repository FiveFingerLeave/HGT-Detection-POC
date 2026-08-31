#!/usr/bin/env python3
"""First-pass mini-chromosome (mChr) candidate classification from the
combined contig metrics table.

This is a provisional, transparent multi-evidence heuristic, not a final
call: per isolate, contigs at/above `--core-min-length-bp` define a "core"
reference (length-weighted mean GC, repeat fraction and gene density).
Smaller contigs are then scored against that reference. Thresholds are
dataset-dependent and must be recalibrated as more isolates and evidence
types (telomeres, synteny, secretome, coverage) are added; see
docs/decisions.md.

Default size thresholds reflect published M. oryzae biology: mini-
chromosomes are supernumerary, repeat-rich, gene-poor accessory
chromosomes roughly a few hundred kb up to ~3 Mb in size. Contigs strictly
between `--small-max-length-bp` and `--core-min-length-bp` fall into an
intentional gray zone ("uncertain") rather than being forced into either
class.

Classes (as defined in Starfish_und_MiniChromosomen_Analyseplan.md, C3):
- core_like: contig length at/above the core threshold.
- mChr_candidate: below the small-contig threshold with >=2 independent
  lines of supporting evidence (repeat-enriched, GC-deviant, gene-poor).
- accessory_candidate: below the small-contig threshold with exactly 1
  supporting evidence line.
- uncertain: below the small-contig threshold with no supporting evidence,
  or in the gray zone between the two length thresholds, or when no core
  reference could be established for the isolate (e.g. a single-contig
  assembly, or all contigs are non-nuclear).
- excluded_non_nuclear: assembly_unit indicates the mitochondrial/
  non-nuclear genome; not part of the core/accessory/mChr scheme.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

COLUMNS = [
    "isolate_id",
    "contig_id",
    "length_bp",
    "gc_fraction",
    "gene_count",
    "repeat_fraction",
    "assembly_unit",
    "role",
    "classification",
]


def _parse_optional_float(value: str) -> float | None:
    return float(value) if value != "" else None


def read_rows(path: Path) -> list[dict[str, object]]:
    with path.open("r", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        rows = []
        for raw in reader:
            rows.append(
                {
                    "isolate_id": raw["isolate_id"],
                    "contig_id": raw["contig_id"],
                    "length_bp": int(raw["length_bp"]),
                    "gc_fraction": _parse_optional_float(raw["gc_fraction"]),
                    "gene_count": _parse_optional_float(raw["gene_count"]),
                    "repeat_fraction": _parse_optional_float(raw["repeat_fraction"]),
                    "assembly_unit": raw["assembly_unit"],
                    "role": raw["role"],
                }
            )
        return rows


def _weighted_mean(values_and_weights: list[tuple[float, float]]) -> float | None:
    total_weight = sum(weight for _, weight in values_and_weights)
    if total_weight == 0:
        return None
    return sum(value * weight for value, weight in values_and_weights) / total_weight


def classify_contigs(
    rows: list[dict[str, object]],
    core_min_length_bp: int,
    small_max_length_bp: int,
    gc_deviation_threshold: float,
    repeat_enrichment_threshold: float,
    gene_density_ratio_threshold: float,
) -> list[dict[str, object]]:
    rows_by_isolate: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        rows_by_isolate.setdefault(row["isolate_id"], []).append(row)

    classified: list[dict[str, object]] = []
    for isolate_rows in rows_by_isolate.values():
        nuclear_rows = [r for r in isolate_rows if r["assembly_unit"] != "non-nuclear"]
        core_rows = [r for r in nuclear_rows if r["length_bp"] >= core_min_length_bp]

        reference_gc = _weighted_mean(
            [(r["gc_fraction"], r["length_bp"]) for r in core_rows if r["gc_fraction"] is not None]
        )
        reference_repeat = _weighted_mean(
            [
                (r["repeat_fraction"], r["length_bp"])
                for r in core_rows
                if r["repeat_fraction"] is not None
            ]
        )
        reference_gene_density = _weighted_mean(
            [
                (r["gene_count"] / r["length_bp"], r["length_bp"])
                for r in core_rows
                if r["gene_count"] is not None
            ]
        )
        has_reference = any(
            ref is not None for ref in (reference_gc, reference_repeat, reference_gene_density)
        )

        for row in isolate_rows:
            row = dict(row)
            if row["assembly_unit"] == "non-nuclear":
                row["classification"] = "excluded_non_nuclear"
            elif row["length_bp"] >= core_min_length_bp:
                row["classification"] = "core_like"
            elif row["length_bp"] > small_max_length_bp or not has_reference:
                row["classification"] = "uncertain"
            else:
                evidence = 0
                if (
                    reference_repeat is not None
                    and row["repeat_fraction"] is not None
                    and row["repeat_fraction"] - reference_repeat >= repeat_enrichment_threshold
                ):
                    evidence += 1
                if (
                    reference_gc is not None
                    and row["gc_fraction"] is not None
                    and abs(row["gc_fraction"] - reference_gc) >= gc_deviation_threshold
                ):
                    evidence += 1
                if (
                    reference_gene_density is not None
                    and reference_gene_density > 0
                    and row["gene_count"] is not None
                ):
                    density = row["gene_count"] / row["length_bp"]
                    if density / reference_gene_density <= gene_density_ratio_threshold:
                        evidence += 1

                if evidence >= 2:
                    row["classification"] = "mChr_candidate"
                elif evidence == 1:
                    row["classification"] = "accessory_candidate"
                else:
                    row["classification"] = "uncertain"

            classified.append(row)

    return classified


def write_rows(rows: list[dict[str, object]], output: Path) -> None:
    with output.open("w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(COLUMNS)
        for row in rows:
            values = []
            for col in COLUMNS:
                value = row[col]
                if isinstance(value, float):
                    values.append(f"{value:.4f}")
                elif value is None:
                    values.append("")
                else:
                    values.append(str(value))
            writer.writerow(values)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Combined contig_metrics TSV")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV with classification")
    parser.add_argument("--core-min-length-bp", type=int, default=4_000_000)
    parser.add_argument("--small-max-length-bp", type=int, default=3_000_000)
    parser.add_argument("--gc-deviation-threshold", type=float, default=0.03)
    parser.add_argument("--repeat-enrichment-threshold", type=float, default=0.10)
    parser.add_argument("--gene-density-ratio-threshold", type=float, default=0.5)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    rows = read_rows(args.input)
    classified = classify_contigs(
        rows,
        args.core_min_length_bp,
        args.small_max_length_bp,
        args.gc_deviation_threshold,
        args.repeat_enrichment_threshold,
        args.gene_density_ratio_threshold,
    )
    write_rows(classified, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
