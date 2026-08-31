from pathlib import Path

import pytest

from normalize_fasta_headers import normalize_fasta_headers


def test_normalize_fasta_headers_writes_prefixed_ids_and_mapping(tmp_path: Path) -> None:
    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(
        ">contig1 description one\nACGTACGT\n>contig2\nTTTTGGGG\nCCCC\n"
    )
    output_fasta = tmp_path / "output.fasta"
    mapping_tsv = tmp_path / "mapping.tsv"

    normalize_fasta_headers(input_fasta, "GCA_TEST", output_fasta, mapping_tsv)

    assert output_fasta.read_text() == (
        ">GCA_TEST_contig1\nACGTACGT\n>GCA_TEST_contig2\nTTTTGGGG\nCCCC\n"
    )
    assert mapping_tsv.read_text() == (
        "old_id\tnew_id\n"
        "contig1\tGCA_TEST_contig1\n"
        "contig2\tGCA_TEST_contig2\n"
    )


def test_normalize_fasta_headers_rejects_empty_header(tmp_path: Path) -> None:
    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(">\nACGTACGT\n")

    with pytest.raises(ValueError, match="empty FASTA header"):
        normalize_fasta_headers(
            input_fasta, "GCA_TEST", tmp_path / "out.fasta", tmp_path / "map.tsv"
        )


def test_normalize_fasta_headers_rejects_duplicate_ids(tmp_path: Path) -> None:
    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(">contig1\nACGT\n>contig1\nTTTT\n")

    with pytest.raises(ValueError, match="duplicate"):
        normalize_fasta_headers(
            input_fasta, "GCA_TEST", tmp_path / "out.fasta", tmp_path / "map.tsv"
        )
