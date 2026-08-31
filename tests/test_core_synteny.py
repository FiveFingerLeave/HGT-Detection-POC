from pathlib import Path

from compute_core_synteny import compute_core_synteny, read_contig_ids


def _paf_line(qname, qlen, qstart, qend, tname, tlen, tstart, tend, mapq=60):
    matches = qend - qstart
    aln_len = matches
    return f"{qname}\t{qlen}\t{qstart}\t{qend}\t+\t{tname}\t{tlen}\t{tstart}\t{tend}\t{matches}\t{aln_len}\t{mapq}\n"


def test_small_contig_with_no_alignment_has_zero_coverage(tmp_path: Path) -> None:
    paf = tmp_path / "self.paf"
    paf.write_text("")

    coverage = compute_core_synteny(paf, ["small1", "core1"], core_min_length_bp=4_000_000)

    assert coverage == {"small1": 0.0, "core1": 0.0}


def test_ignores_self_hits(tmp_path: Path) -> None:
    paf = tmp_path / "self.paf"
    paf.write_text(_paf_line("core1", 5_000_000, 0, 5_000_000, "core1", 5_000_000, 0, 5_000_000))

    coverage = compute_core_synteny(paf, ["core1"], core_min_length_bp=4_000_000)

    assert coverage["core1"] == 0.0


def test_ignores_alignment_to_non_core_sized_target(tmp_path: Path) -> None:
    paf = tmp_path / "self.paf"
    paf.write_text(
        _paf_line("small1", 100_000, 0, 100_000, "small2", 100_000, 0, 100_000)
    )

    coverage = compute_core_synteny(
        paf, ["small1", "small2"], core_min_length_bp=4_000_000
    )

    assert coverage["small1"] == 0.0


def test_full_length_alignment_to_core_contig_gives_full_coverage(tmp_path: Path) -> None:
    paf = tmp_path / "self.paf"
    paf.write_text(
        _paf_line("small1", 100_000, 0, 100_000, "core1", 5_000_000, 1_000_000, 1_100_000)
    )

    coverage = compute_core_synteny(paf, ["small1", "core1"], core_min_length_bp=4_000_000)

    assert coverage["small1"] == 1.0


def test_overlapping_alignment_blocks_are_merged_not_double_counted(tmp_path: Path) -> None:
    paf = tmp_path / "self.paf"
    lines = _paf_line(
        "small1", 100_000, 0, 60_000, "core1", 5_000_000, 0, 60_000
    ) + _paf_line("small1", 100_000, 40_000, 100_000, "core1", 5_000_000, 200_000, 260_000)
    paf.write_text(lines)

    coverage = compute_core_synteny(paf, ["small1", "core1"], core_min_length_bp=4_000_000)

    # merged union of [0,60000) and [40000,100000) = [0,100000) = full length
    assert coverage["small1"] == 1.0


def test_read_contig_ids_preserves_fasta_order(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">contig2 desc\nACGT\n>contig1\nACGT\n")

    assert read_contig_ids(fasta) == ["contig2", "contig1"]
