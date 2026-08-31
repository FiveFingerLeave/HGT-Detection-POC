from pathlib import Path

from validate_fasta_gff_ids import validate_fasta_gff_ids


def test_validate_fasta_gff_ids_returns_empty_when_compatible(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">contig1 description\nACGT\n>contig2\nTTTT\n")
    gff = tmp_path / "genome.gff"
    gff.write_text(
        "##gff-version 3\n"
        "##sequence-region contig1 1 4\n"
        "contig1\tGenbank\tregion\t1\t4\t.\t+\t.\tID=contig1\n"
        "contig2\tGenbank\tregion\t1\t4\t.\t+\t.\tID=contig2\n"
    )

    assert validate_fasta_gff_ids(fasta, gff) == []


def test_validate_fasta_gff_ids_reports_missing_seqids(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">contig1\nACGT\n")
    gff = tmp_path / "genome.gff"
    gff.write_text(
        "##sequence-region contig2 1 4\n"
        "contig2\tGenbank\tregion\t1\t4\t.\t+\t.\tID=contig2\n"
    )

    assert validate_fasta_gff_ids(fasta, gff) == ["contig2"]
