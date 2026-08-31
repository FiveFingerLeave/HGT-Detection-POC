from pathlib import Path

from compute_reference_synteny import compute_reference_synteny


def _paf_line(qname, qstart, qend, tname, tstart, tend, mapq=60):
    matches = qend - qstart
    return f"{qname}\t999999\t{qstart}\t{qend}\t+\t{tname}\t999999\t{tstart}\t{tend}\t{matches}\t{matches}\t{mapq}\n"


def _write_core_contigs(path: Path, contig_ids: list[str]) -> Path:
    path.write_text("\n".join(contig_ids) + "\n")
    return path


def test_contig_matching_enough_references_is_consensus(tmp_path: Path) -> None:
    query_fasta = tmp_path / "query.fasta"
    query_fasta.write_text(f">small\n{'A' * 100_000}\n")

    paf1 = tmp_path / "ref1.paf"
    paf1.write_text(_paf_line("small", 0, 100_000, "ref1_chr1", 0, 100_000))
    core1 = _write_core_contigs(tmp_path / "ref1_core.txt", ["ref1_chr1"])

    paf2 = tmp_path / "ref2.paf"
    paf2.write_text(_paf_line("small", 0, 100_000, "ref2_chr1", 0, 100_000))
    core2 = _write_core_contigs(tmp_path / "ref2_core.txt", ["ref2_chr1"])

    result = compute_reference_synteny(
        query_fasta,
        reference_pafs=[paf1, paf2],
        reference_core_contig_files=[core1, core2],
        min_coverage_per_reference=0.3,
        consensus_min_hits=2,
    )

    assert result["small"]["reference_core_hits"] == 2
    assert result["small"]["reference_core_consensus"] is True


def test_contig_below_consensus_threshold_is_not_consensus(tmp_path: Path) -> None:
    query_fasta = tmp_path / "query.fasta"
    query_fasta.write_text(f">small\n{'A' * 100_000}\n")

    paf1 = tmp_path / "ref1.paf"
    paf1.write_text(_paf_line("small", 0, 100_000, "ref1_chr1", 0, 100_000))
    core1 = _write_core_contigs(tmp_path / "ref1_core.txt", ["ref1_chr1"])

    paf2 = tmp_path / "ref2.paf"
    paf2.write_text("")  # no alignment at all to reference 2
    core2 = _write_core_contigs(tmp_path / "ref2_core.txt", ["ref2_chr1"])

    result = compute_reference_synteny(
        query_fasta,
        reference_pafs=[paf1, paf2],
        reference_core_contig_files=[core1, core2],
        min_coverage_per_reference=0.3,
        consensus_min_hits=2,
    )

    assert result["small"]["reference_core_hits"] == 1
    assert result["small"]["reference_core_consensus"] is False


def test_alignment_to_non_core_reference_contig_is_ignored(tmp_path: Path) -> None:
    query_fasta = tmp_path / "query.fasta"
    query_fasta.write_text(f">small\n{'A' * 100_000}\n")

    paf1 = tmp_path / "ref1.paf"
    # aligns fully, but to a contig NOT in ref1's core list (e.g. its own mChr)
    paf1.write_text(_paf_line("small", 0, 100_000, "ref1_mchr", 0, 100_000))
    core1 = _write_core_contigs(tmp_path / "ref1_core.txt", ["ref1_chr1"])

    result = compute_reference_synteny(
        query_fasta,
        reference_pafs=[paf1],
        reference_core_contig_files=[core1],
        min_coverage_per_reference=0.3,
        consensus_min_hits=1,
    )

    assert result["small"]["reference_core_hits"] == 0


def test_coverage_below_per_reference_threshold_does_not_count_as_hit(tmp_path: Path) -> None:
    query_fasta = tmp_path / "query.fasta"
    query_fasta.write_text(f">small\n{'A' * 100_000}\n")

    paf1 = tmp_path / "ref1.paf"
    paf1.write_text(_paf_line("small", 0, 10_000, "ref1_chr1", 0, 10_000))  # only 10% covered
    core1 = _write_core_contigs(tmp_path / "ref1_core.txt", ["ref1_chr1"])

    result = compute_reference_synteny(
        query_fasta,
        reference_pafs=[paf1],
        reference_core_contig_files=[core1],
        min_coverage_per_reference=0.3,
        consensus_min_hits=1,
    )

    assert result["small"]["reference_core_hits"] == 0


def test_scattered_hits_across_many_targets_do_not_sum_to_a_hit(tmp_path: Path) -> None:
    # Regression test: a repeat/TE family present in many different core
    # chromosomes can produce small alignment blocks against several
    # different targets that would sum to a large total, even though no
    # single target shows real synteny. Coverage must be taken as the best
    # SINGLE-target coverage, not summed across targets.
    query_fasta = tmp_path / "query.fasta"
    query_fasta.write_text(f">small\n{'A' * 100_000}\n")

    paf1 = tmp_path / "ref1.paf"
    paf1.write_text(
        _paf_line("small", 0, 15_000, "ref1_chr1", 0, 15_000)
        + _paf_line("small", 20_000, 35_000, "ref1_chr2", 0, 15_000)
        + _paf_line("small", 40_000, 55_000, "ref1_chr3", 0, 15_000)
    )
    core1 = _write_core_contigs(
        tmp_path / "ref1_core.txt", ["ref1_chr1", "ref1_chr2", "ref1_chr3"]
    )

    result = compute_reference_synteny(
        query_fasta,
        reference_pafs=[paf1],
        reference_core_contig_files=[core1],
        min_coverage_per_reference=0.3,
        consensus_min_hits=1,
    )

    # summed coverage would be 45% (>= 0.3); best single-target is 15%
    assert result["small"]["reference_core_hits"] == 0


def test_low_mapping_quality_alignment_is_excluded(tmp_path: Path) -> None:
    query_fasta = tmp_path / "query.fasta"
    query_fasta.write_text(f">small\n{'A' * 100_000}\n")

    paf1 = tmp_path / "ref1.paf"
    paf1.write_text(_paf_line("small", 0, 100_000, "ref1_chr1", 0, 100_000, mapq=5))
    core1 = _write_core_contigs(tmp_path / "ref1_core.txt", ["ref1_chr1"])

    result = compute_reference_synteny(
        query_fasta,
        reference_pafs=[paf1],
        reference_core_contig_files=[core1],
        min_coverage_per_reference=0.3,
        consensus_min_hits=1,
        min_mapq=30,
    )

    assert result["small"]["reference_core_hits"] == 0


def test_contig_with_no_alignments_anywhere_has_zero_hits(tmp_path: Path) -> None:
    query_fasta = tmp_path / "query.fasta"
    query_fasta.write_text(f">lonely\n{'A' * 100_000}\n")

    paf1 = tmp_path / "ref1.paf"
    paf1.write_text("")
    core1 = _write_core_contigs(tmp_path / "ref1_core.txt", ["ref1_chr1"])

    result = compute_reference_synteny(
        query_fasta,
        reference_pafs=[paf1],
        reference_core_contig_files=[core1],
        min_coverage_per_reference=0.3,
        consensus_min_hits=1,
    )

    assert result["lonely"]["reference_core_hits"] == 0
    assert result["lonely"]["reference_core_consensus"] is False
