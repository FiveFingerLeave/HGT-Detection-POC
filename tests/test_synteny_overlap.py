from synteny_overlap import (
    aligned_fraction,
    classify_synteny,
    evaluate_regions,
    merge_intervals,
)


def test_merge_intervals_joins_overlapping():
    assert merge_intervals([(0, 100), (50, 150), (200, 300)]) == [[0, 150], [200, 300]]


def test_merge_intervals_empty():
    assert merge_intervals([]) == []


def test_aligned_fraction_full_coverage():
    assert aligned_fraction([[0, 100]], 0, 100) == 1.0


def test_aligned_fraction_partial_coverage():
    assert aligned_fraction([[0, 50]], 0, 100) == 0.5


def test_aligned_fraction_no_coverage():
    assert aligned_fraction([[200, 300]], 0, 100) == 0.0


def test_aligned_fraction_zero_length_region():
    assert aligned_fraction([[0, 100]], 50, 50) == 0.0


def test_classify_syntenic_present_when_region_aligns():
    assert classify_synteny(0.95, 0.9, 0.9) == "syntenic_present"


def test_classify_syntenic_empty_site_when_flanks_align_but_not_region():
    # The mobile-element insertion-polymorphism signature: flanks conserved,
    # cargo (the region itself) absent in the other genome.
    assert classify_synteny(0.0, 0.9, 0.85) == "syntenic_empty_site"


def test_classify_non_syntenic_when_nothing_aligns():
    assert classify_synteny(0.05, 0.1, 0.05) == "non_syntenic_region"


def test_classify_ambiguous_otherwise():
    assert classify_synteny(0.4, 0.6, 0.2) == "ambiguous"


def test_evaluate_regions_end_to_end():
    regions = [
        {"chrom": "chr1", "start": 10_000, "end": 12_000, "region_id": "r1", "region_class": "candidate"}
    ]
    # flanks align (8000-10000 and 12000-14000), the region itself does not.
    paf_intervals = {"chr1": [(7_000, 10_000), (12_000, 15_000)]}

    rows = evaluate_regions(regions, paf_intervals, flank_bp=2000, genome="A", other_genome="B")

    assert len(rows) == 1
    assert rows[0]["classification"] == "syntenic_empty_site"
    assert rows[0]["genome"] == "A"
    assert rows[0]["other_genome"] == "B"


def test_evaluate_regions_filters_by_valid_contigs_upstream():
    # (read_regions filtering is exercised via the CLI/main path; here we
    # just confirm evaluate_regions handles a chrom absent from the PAF
    # gracefully, i.e. reports zero alignment rather than crashing.)
    regions = [
        {"chrom": "chr_unmapped", "start": 0, "end": 100, "region_id": "r2", "region_class": "core"}
    ]
    rows = evaluate_regions(regions, {}, flank_bp=500, genome="A", other_genome="B")

    assert rows[0]["region_aligned_fraction"] == 0.0
    assert rows[0]["classification"] == "non_syntenic_region"
