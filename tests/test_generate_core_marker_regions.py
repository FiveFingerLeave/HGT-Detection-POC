from generate_core_marker_regions import generate_core_windows


def test_generate_core_windows_uses_largest_contig():
    contig_lengths = {"small": 10_000, "big": 1_000_000}

    regions = generate_core_windows(contig_lengths, {}, "iso1", window_size=10_000, n_windows=3)

    assert all(r[0] == "big" for r in regions)
    assert len(regions) == 3


def test_generate_core_windows_skips_excluded_overlap():
    contig_lengths = {"big": 1_000_000}
    stride = 1_000_000 // (3 + 1)  # matches internal stride for n_windows=3
    excluded = {"big": [(stride, stride + 10_000)]}  # covers the first window exactly

    regions = generate_core_windows(contig_lengths, excluded, "iso1", window_size=10_000, n_windows=3)

    starts = [r[1] for r in regions]
    assert stride not in starts


def test_generate_core_windows_empty_when_no_contigs():
    assert generate_core_windows({}, {}, "iso1", 10_000, 3) == []
