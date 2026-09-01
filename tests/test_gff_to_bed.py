from pathlib import Path

from gff_to_bed import gff_to_regions


def test_gff_to_regions_extracts_mrna_only(tmp_path: Path) -> None:
    gff = tmp_path / "genome_YR.filt.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tMetaEuk\tgene\t100\t200\t.\t+\t.\tID=gene1\n"
        "chr1\tMetaEuk\tmRNA\t100\t200\t.\t+\t.\tID=mrna1;Parent=gene1\n"
        "chr1\tMetaEuk\texon\t100\t150\t.\t+\t.\tID=exon1\n"
        "chr2\tMetaEuk\tmRNA\t500\t600\t.\t+\t.\tID=mrna2\n"
    )

    regions = gff_to_regions(gff, "iso1", "candidate", "starship")

    assert regions == [
        ("chr1", 99, 200, "starship_iso1_1", "candidate"),
        ("chr2", 499, 600, "starship_iso1_2", "candidate"),
    ]
