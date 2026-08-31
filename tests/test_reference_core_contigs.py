from pathlib import Path

from list_reference_core_contigs import list_core_contigs


def test_selects_contigs_at_or_above_length_threshold(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">big\n{'A' * 5_000_000}\n>small\n{'A' * 100_000}\n")

    core = list_core_contigs(fasta, core_min_length_bp=4_000_000)

    assert core == ["big"]


def test_excludes_non_nuclear_contig_even_if_large(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">chr1\n{'A' * 5_000_000}\n>mt\n{'A' * 4_500_000}\n")
    sequence_report = tmp_path / "sequence_report.jsonl"
    sequence_report.write_text(
        '{"genbankAccession": "chr1", "assemblyUnit": "Primary Assembly"}\n'
        '{"genbankAccession": "mt", "assemblyUnit": "non-nuclear"}\n'
    )

    core = list_core_contigs(
        fasta, core_min_length_bp=4_000_000, sequence_report=sequence_report
    )

    assert core == ["chr1"]


def test_excludes_explicitly_listed_accessory_contig(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">chr1\n{'A' * 5_000_000}\n>mchr\n{'A' * 4_200_000}\n")

    core = list_core_contigs(
        fasta, core_min_length_bp=4_000_000, exclude={"mchr"}
    )

    assert core == ["chr1"]


def test_works_when_all_contigs_report_unplaced_scaffold_role(tmp_path: Path) -> None:
    # Contig-level assemblies (e.g. CD156, US71) mark every contig as
    # 'unplaced-scaffold' regardless of size - length must still decide.
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">contig1\n{'A' * 6_000_000}\n>contig2\n{'A' * 50_000}\n")
    sequence_report = tmp_path / "sequence_report.jsonl"
    sequence_report.write_text(
        '{"genbankAccession": "contig1", "assemblyUnit": "Primary Assembly", "role": "unplaced-scaffold"}\n'
        '{"genbankAccession": "contig2", "assemblyUnit": "Primary Assembly", "role": "unplaced-scaffold"}\n'
    )

    core = list_core_contigs(
        fasta, core_min_length_bp=4_000_000, sequence_report=sequence_report
    )

    assert core == ["contig1"]
