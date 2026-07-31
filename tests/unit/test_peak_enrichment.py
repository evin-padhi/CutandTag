from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

from peak_enrichment import Interval, count_overlapping_foreground, parse_intervals


def test_parse_intervals_rejects_unknown_chromosome_and_invalid_coordinates(tmp_path):
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t20\t10\nchrX\t0\t5\n", encoding="utf-8")

    with pytest.raises(ValueError, match="start must be less than end"):
        parse_intervals(peaks, {"chr1": 100})


def test_parse_intervals_rejects_unknown_chromosome(tmp_path):
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chrX\t0\t5\n", encoding="utf-8")

    with pytest.raises(ValueError, match="unknown chromosome"):
        parse_intervals(peaks, {"chr1": 100})


def test_overlap_counts_each_foreground_interval_once():
    foreground = [Interval("chr1", 10, 20), Interval("chr1", 30, 40)]
    reference = [Interval("chr1", 15, 17), Interval("chr1", 16, 18)]

    assert count_overlapping_foreground(foreground, reference) == 1


def test_gc_fraction_reports_fraction_of_gc_bases():
    from peak_enrichment import gc_fraction

    fasta = {"chr1": "AACCGGTT"}
    interval = Interval("chr1", 2, 6)

    assert gc_fraction(interval, fasta) == 1.0


def test_length_gc_background_preserves_width_chromosome_and_gc():
    import random

    from peak_enrichment import gc_fraction, sample_background

    fasta = {"chr1": "ACGT" * 100}
    foreground = [Interval("chr1", 8, 20)]
    background = sample_background(
        foreground,
        {"chr1": 400},
        fasta,
        [],
        "length_gc_matched",
        random.Random(7),
        max_attempts=1000,
    )

    assert background is not None
    assert background[0].chrom == "chr1"
    assert background[0].end - background[0].start == 12
    assert abs(gc_fraction(background[0], fasta) - gc_fraction(foreground[0], fasta)) <= 0.02


def test_blacklisted_candidate_is_not_sampled():
    import random

    from peak_enrichment import sample_background

    foreground = [Interval("chr1", 10, 20)]
    blacklist = [Interval("chr1", 30, 40)]
    background = sample_background(
        foreground,
        {"chr1": 100},
        {"chr1": "A" * 100},
        blacklist,
        "length_matched",
        random.Random(3),
        max_attempts=1000,
    )

    assert background is not None
    assert background[0] != Interval("chr1", 30, 40)


def test_calculate_enrichment_reports_ratio_and_upper_tail_p_value():
    from peak_enrichment import calculate_enrichment

    foreground = [Interval("chr1", 10, 20), Interval("chr1", 40, 50)]
    reference = [Interval("chr1", 10, 20)]
    rows = calculate_enrichment(
        foreground,
        reference,
        {"chr1": 100},
        {"chr1": "ACGT" * 25},
        [],
        permutations=20,
        seed=11,
        gc_tolerance=0.02,
    )

    assert {row["background_model"] for row in rows} == {
        "random",
        "length_matched",
        "gc_matched",
        "length_gc_matched",
    }
    assert all(row["observed_overlap_count"] == 1 for row in rows)
    assert all(0 <= row["empirical_p_value"] <= 1 for row in rows)
    assert all(row["permutations_succeeded"] <= 20 for row in rows)


def test_calculate_enrichment_emits_status_rows_for_empty_inputs():
    from peak_enrichment import calculate_enrichment

    no_foreground_rows = calculate_enrichment(
        [],
        [Interval("chr1", 10, 20)],
        {"chr1": 100},
        {"chr1": "ACGT" * 25},
        [],
        permutations=5,
        seed=7,
        gc_tolerance=0.02,
    )
    no_reference_rows = calculate_enrichment(
        [Interval("chr1", 10, 20)],
        [],
        {"chr1": 100},
        {"chr1": "ACGT" * 25},
        [],
        permutations=5,
        seed=7,
        gc_tolerance=0.02,
    )

    assert {row["status"] for row in no_foreground_rows} == {"no_foreground_peaks"}
    assert {row["status"] for row in no_reference_rows} == {"no_reference_peaks"}
    assert all(row["observed_overlap_count"] is None for row in no_foreground_rows)
    assert all(row["observed_overlap_count"] is None for row in no_reference_rows)
