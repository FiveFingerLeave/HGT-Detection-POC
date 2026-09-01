#!/usr/bin/env python3
"""Convert one or more GFF files (each labelled with the genome/reference
they came from) into a single BED file of candidate regions, for
candidate_regions.bed (Section 1/3 of the POC guideline).

Only 'mRNA' features are taken (one row per predicted transcript, matching
starfish's own filt.gff convention) to avoid duplicate/nested coordinates
from gene+mRNA+exon+CDS all describing the same locus.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def gff_to_regions(
    gff: Path, genome_id: str, region_class: str, id_prefix: str, feature_type: str = "mRNA"
) -> list[tuple[str, int, int, str, str]]:
    regions = []
    n = 0
    with gff.open("r") as infile:
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != feature_type:
                continue
            chrom, start, end = fields[0], int(fields[3]) - 1, int(fields[4])
            n += 1
            region_id = f"{id_prefix}_{genome_id}_{n}"
            regions.append((chrom, start, end, region_id, region_class))
    return regions


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gff",
        action="append",
        required=True,
        metavar="GENOME_ID=PATH",
        dest="gff_pairs",
        help="Repeatable: genome_id=path/to/genome_YR.filt.gff",
    )
    parser.add_argument("--region-class", default="candidate")
    parser.add_argument("--id-prefix", default="region")
    parser.add_argument("--feature-type", default="mRNA")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    all_regions = []
    for pair in args.gff_pairs:
        genome_id, path = pair.split("=", 1)
        all_regions.extend(
            gff_to_regions(Path(path), genome_id, args.region_class, args.id_prefix, args.feature_type)
        )

    with args.output.open("w") as outfile:
        for chrom, start, end, region_id, region_class in all_regions:
            outfile.write(f"{chrom}\t{start}\t{end}\t{region_id}\t{region_class}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
