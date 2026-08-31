from pathlib import Path

import pytest

from normalize_gff_seqids import normalize_gff_seqids


def _write_mapping(tmp_path: Path) -> Path:
    mapping_tsv = tmp_path / "mapping.tsv"
    mapping_tsv.write_text("old_id\tnew_id\ncontig1\tGCA_TEST_contig1\n")
    return mapping_tsv


def test_normalize_gff_seqids_rewrites_pragma_and_feature_lines(tmp_path: Path) -> None:
    input_gff = tmp_path / "input.gff"
    input_gff.write_text(
        "##gff-version 3\n"
        "##sequence-region contig1 1 1000\n"
        "contig1\tGenbank\tregion\t1\t1000\t.\t+\t.\tID=contig1:1..1000\n"
    )
    mapping_tsv = _write_mapping(tmp_path)
    output_gff = tmp_path / "output.gff"

    normalize_gff_seqids(input_gff, mapping_tsv, output_gff)

    assert output_gff.read_text() == (
        "##gff-version 3\n"
        "##sequence-region GCA_TEST_contig1 1 1000\n"
        "GCA_TEST_contig1\tGenbank\tregion\t1\t1000\t.\t+\t.\tID=contig1:1..1000\n"
    )


def test_normalize_gff_seqids_rejects_seqid_without_mapping(tmp_path: Path) -> None:
    input_gff = tmp_path / "input.gff"
    input_gff.write_text("contig2\tGenbank\tregion\t1\t1000\t.\t+\t.\tID=contig2:1..1000\n")
    mapping_tsv = _write_mapping(tmp_path)

    with pytest.raises(ValueError, match="no entry in the id mapping"):
        normalize_gff_seqids(input_gff, mapping_tsv, tmp_path / "output.gff")
