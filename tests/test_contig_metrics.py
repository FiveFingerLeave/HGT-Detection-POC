from pathlib import Path

from calculate_contig_metrics import calculate_contig_metrics


def test_calculate_contig_metrics_without_gff(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">contig1\nGGCC\nAATT\n>contig2\nNNNNGGGG\n")

    rows = calculate_contig_metrics(fasta, "ISO1")

    assert rows == [
        {
            "isolate_id": "ISO1",
            "contig_id": "contig1",
            "length_bp": 8,
            "gc_fraction": 0.5,
            "gene_count": "",
            "repeat_fraction": "",
            "assembly_unit": "",
            "role": "",
        },
        {
            "isolate_id": "ISO1",
            "contig_id": "contig2",
            "length_bp": 8,
            "gc_fraction": 1.0,
            "gene_count": "",
            "repeat_fraction": "",
            "assembly_unit": "",
            "role": "",
        },
    ]


def test_calculate_contig_metrics_with_gff_adds_gene_count(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">contig1\nACGT\n>contig2\nACGT\n")
    gff = tmp_path / "genome.gff"
    gff.write_text(
        "##gff-version 3\n"
        "contig1\tGenbank\tgene\t1\t4\t.\t+\t.\tID=gene1\n"
        "contig1\tGenbank\tmRNA\t1\t4\t.\t+\t.\tID=mrna1;Parent=gene1\n"
        "contig1\tGenbank\tgene\t1\t4\t.\t+\t.\tID=gene2\n"
    )

    rows = calculate_contig_metrics(fasta, "ISO1", gff)

    assert rows[0]["gene_count"] == 2
    assert rows[1]["gene_count"] == 0


def test_calculate_contig_metrics_all_n_contig_has_empty_gc(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">contig1\nNNNN\n")

    rows = calculate_contig_metrics(fasta, "ISO1")

    assert rows[0]["gc_fraction"] == ""


def test_calculate_contig_metrics_flags_non_nuclear_contig(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">ISO1-CP1\nACGT\n>ISO1-MT1\nACGT\n")
    id_map = tmp_path / "id_map.tsv"
    id_map.write_text("old_id\tnew_id\nCP1\tISO1-CP1\nMT1\tISO1-MT1\n")
    sequence_report = tmp_path / "sequence_report.jsonl"
    sequence_report.write_text(
        '{"genbankAccession": "CP1", "assemblyUnit": "Primary Assembly", "role": "assembled-molecule"}\n'
        '{"genbankAccession": "MT1", "assemblyUnit": "non-nuclear", "role": "assembled-molecule"}\n'
    )

    rows = calculate_contig_metrics(
        fasta, "ISO1", sequence_report=sequence_report, id_map=id_map
    )

    assert rows[0]["assembly_unit"] == "Primary Assembly"
    assert rows[1]["assembly_unit"] == "non-nuclear"


def test_calculate_contig_metrics_flags_unplaced_scaffold_role(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">ISO1-A\nACGT\n>ISO1-B\nACGT\n")
    id_map = tmp_path / "id_map.tsv"
    id_map.write_text("old_id\tnew_id\nA\tISO1-A\nB\tISO1-B\n")
    sequence_report = tmp_path / "sequence_report.jsonl"
    sequence_report.write_text(
        '{"genbankAccession": "A", "assemblyUnit": "Primary Assembly", "role": "assembled-molecule"}\n'
        '{"genbankAccession": "B", "assemblyUnit": "Primary Assembly", "role": "unplaced-scaffold"}\n'
    )

    rows = calculate_contig_metrics(
        fasta, "ISO1", sequence_report=sequence_report, id_map=id_map
    )

    assert rows[0]["role"] == "assembled-molecule"
    assert rows[1]["role"] == "unplaced-scaffold"


def test_calculate_contig_metrics_adds_repeat_fraction_from_masked_fasta(
    tmp_path: Path,
) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">contig1\nACGTACGT\n>contig2\nACGTACGT\n")
    masked_fasta = tmp_path / "genome.masked.fasta"
    masked_fasta.write_text(">contig1\nACGTacgt\n>contig2\nACGTACGT\n")

    rows = calculate_contig_metrics(fasta, "ISO1", masked_fasta=masked_fasta)

    assert rows[0]["repeat_fraction"] == 0.5
    assert rows[1]["repeat_fraction"] == 0.0
