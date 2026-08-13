import csv
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

from peak_enrichment import Interval, count_overlapping_foreground, parse_intervals


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_enrichment_fixture(tmp_path: Path) -> Path:
    from peak_enrichment import main

    outdir = tmp_path / "enrichment"
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + ("ACGT" * 50) + "\n", encoding="utf-8")

    foreground_dir = tmp_path / "foreground"
    foreground_dir.mkdir()
    (foreground_dir / "ctcf.bed").write_text(
        "chr1\t10\t20\nchr1\t40\t50\n",
        encoding="utf-8",
    )

    reference_dir = tmp_path / "references"
    reference_dir.mkdir()
    (reference_dir / "ctcf_called.bed").write_text("chr1\t10\t20\n", encoding="utf-8")
    (reference_dir / "ctcf_prior.bed").write_text("chr1\t12\t18\nchr1\t70\t80\n", encoding="utf-8")

    foreground_manifest = tmp_path / "foregrounds.tsv"
    foreground_manifest.write_text(
        "foreground_id\tforeground_tf\tpeak_file\n"
        "fg_ctcf\tCTCF\tforeground/ctcf.bed\n",
        encoding="utf-8",
    )
    reference_manifest = tmp_path / "references.tsv"
    reference_manifest.write_text(
        "reference_id\ttf\treference_type\tpeak_file\n"
        "fg_ctcf\tCTCF\tcalled_tf\treferences/ctcf_called.bed\n"
        "ref_ctcf_prior\tCTCF\tchipseq\treferences/ctcf_prior.bed\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "--foreground-manifest",
            str(foreground_manifest),
            "--reference-manifest",
            str(reference_manifest),
            "--fasta",
            str(fasta),
            "--outdir",
            str(outdir),
            "--permutations",
            "8",
            "--seed",
            "1729",
            "--gc-tolerance",
            "0.02",
        ]
    )

    assert exit_code == 0
    return outdir


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


def test_parse_intervals_reads_gzip_bed(tmp_path):
    import gzip

    peaks = tmp_path / "peaks.bed.gz"
    with gzip.open(peaks, "wt", encoding="utf-8") as handle:
        handle.write("chr1\t10\t20\tpeak-1\n")

    assert parse_intervals(peaks, {"chr1": 100}) == [Interval("chr1", 10, 20)]


def test_load_fasta_chrom_sizes_streams_wrapped_records(tmp_path, monkeypatch):
    from peak_enrichment import load_fasta_chrom_sizes

    fasta = tmp_path / "reference.fa"
    fasta.write_text(
        ">chr1 first record\nACG\nTT\n\n>chr2\nA\nCCC\n",
        encoding="utf-8",
    )

    def reject_full_file_read(*args, **kwargs):
        raise AssertionError("chromosome-size loading must stream the FASTA")

    monkeypatch.setattr(Path, "read_text", reject_full_file_read)

    assert load_fasta_chrom_sizes(fasta) == {"chr1": 5, "chr2": 4}


def test_load_fasta_chrom_sizes_rejects_duplicate_record_names(tmp_path):
    from peak_enrichment import load_fasta_chrom_sizes

    fasta = tmp_path / "duplicate.fa"
    fasta.write_text(
        ">chr1 first\nACGT\n>chr1 second\nTGCA\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate FASTA record name 'chr1'"):
        load_fasta_chrom_sizes(fasta)


@pytest.mark.parametrize("header", [">", ">   "])
def test_load_fasta_chrom_sizes_rejects_empty_record_names(tmp_path, header):
    from peak_enrichment import load_fasta_chrom_sizes

    fasta = tmp_path / "empty-name.fa"
    fasta.write_text(f"{header}\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="missing a sequence name"):
        load_fasta_chrom_sizes(fasta)


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


def test_gc_tolerance_changes_candidate_acceptance():
    from peak_enrichment import sample_background

    class SequenceRNG:
        def __init__(self, starts):
            self.starts = list(starts)

        def choice(self, values):
            return values[0]

        def randrange(self, start, stop=None, step=1):
            if stop is None:
                stop = start
                start = 0
            return self.starts.pop(0)

    foreground = [Interval("chr1", 1, 5)]
    fasta = {"chr1": "AAAACCCC"}
    chrom_sizes = {"chr1": 8}

    strict_background = sample_background(
        foreground,
        chrom_sizes,
        fasta,
        [],
        "length_gc_matched",
        SequenceRNG([3, 0]),
        max_attempts=2,
        gc_tolerance=0.01,
    )
    loose_background = sample_background(
        foreground,
        chrom_sizes,
        fasta,
        [],
        "length_gc_matched",
        SequenceRNG([3, 0]),
        max_attempts=2,
        gc_tolerance=0.30,
    )

    assert strict_background is None
    assert loose_background == [Interval("chr1", 0, 4)]


def test_gc_matched_preserves_each_foreground_chromosome_while_sampling_widths():
    import random

    from peak_enrichment import sample_background

    class LastChoiceRNG(random.Random):
        def __init__(self):
            super().__init__(0)
            self.starts = iter((0, 30))

        def choice(self, values):
            return values[-1]

        def randrange(self, start, stop=None, step=1):
            return next(self.starts)

    foreground = [Interval("chr1", 5, 15), Interval("chr2", 5, 25)]
    background = sample_background(
        foreground,
        {"chr1": 100, "chr2": 100},
        {"chr1": "A" * 100, "chr2": "A" * 100},
        [],
        "gc_matched",
        LastChoiceRNG(),
        max_attempts=4,
    )

    assert background is not None
    assert [interval.chrom for interval in background] == ["chr1", "chr2"]
    assert [interval.width for interval in background] == [20, 20]


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


def test_calculate_enrichment_marks_partial_permutations_insufficient(monkeypatch):
    import peak_enrichment

    calls_by_model = {}

    def sample_once_then_fail(
        foreground,
        chrom_sizes,
        fasta_sequences,
        blacklist,
        model,
        rng,
        max_attempts,
        gc_tolerance,
    ):
        calls_by_model[model] = calls_by_model.get(model, 0) + 1
        if calls_by_model[model] == 1:
            return [Interval("chr1", 50, 60)]
        return None

    monkeypatch.setattr(peak_enrichment, "sample_background", sample_once_then_fail)

    rows = peak_enrichment.calculate_enrichment(
        [Interval("chr1", 10, 20)],
        [Interval("chr1", 10, 20)],
        {"chr1": 100},
        {"chr1": "A" * 100},
        [],
        permutations=3,
        seed=11,
        gc_tolerance=0.02,
    )

    assert {row["status"] for row in rows} == {"insufficient_background"}
    assert {row["permutations_succeeded"] for row in rows} == {1}
    assert all(row["observed_overlap_count"] == 1 for row in rows)
    for field in (
        "null_mean_overlap",
        "null_sd_overlap",
        "enrichment_ratio",
        "empirical_p_value",
    ):
        assert all(row[field] is None for row in rows)


def test_calculate_enrichment_marks_zero_null_mean_non_success(monkeypatch):
    import peak_enrichment

    monkeypatch.setattr(
        peak_enrichment,
        "sample_background",
        lambda *args, **kwargs: [Interval("chr1", 50, 60)],
    )

    rows = peak_enrichment.calculate_enrichment(
        [Interval("chr1", 10, 20)],
        [Interval("chr1", 10, 20)],
        {"chr1": 100},
        {"chr1": "A" * 100},
        [],
        permutations=3,
        seed=11,
        gc_tolerance=0.02,
    )

    assert {row["status"] for row in rows} == {"zero_null_mean"}
    assert {row["null_mean_overlap"] for row in rows} == {0}
    assert {row["null_sd_overlap"] for row in rows} == {0.0}
    assert all(row["enrichment_ratio"] is None for row in rows)
    assert all(row["empirical_p_value"] is None for row in rows)


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


def test_cli_resolves_chipseq_peak_paths_relative_to_manifest(tmp_path):
    from peak_enrichment import load_reference_manifest

    manifest = tmp_path / "chipseq.csv"
    peaks = tmp_path / "prior" / "ctcf.bed"
    peaks.parent.mkdir()
    peaks.write_text("chr1\t10\t20\n", encoding="utf-8")
    manifest.write_text(
        "reference_id,tf,peak_file\nctcf,CTCF,prior/ctcf.bed\n",
        encoding="utf-8",
    )

    references = load_reference_manifest(manifest, {"chr1": 100})

    assert references[0].peak_file == peaks.resolve()
    assert references[0].reference_type == "chipseq"


def test_public_chipseq_validation_parses_csv_and_forces_chipseq_type(
    tmp_path,
    capsys,
):
    import json

    from peak_enrichment import main

    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + ("ACGT" * 25) + "\n", encoding="utf-8")
    peaks = tmp_path / "prior,quoted.bed"
    peaks.write_text("chr1\t10\t20\n", encoding="utf-8")
    manifest = tmp_path / "chipseq.csv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["reference_id", "tf", "peak_file", "reference_type"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "reference_id": "public_ctcf",
                "tf": "CTCF, public reference",
                "peak_file": peaks.name,
                "reference_type": "called_tf",
            }
        )

    exit_code = main(
        [
            "validate-chipseq-manifest",
            "--manifest",
            str(manifest),
            "--fasta",
            str(fasta),
        ]
    )

    assert exit_code == 0
    rows = json.loads(capsys.readouterr().out)
    assert rows == [
        {
            "reference_id": "public_ctcf",
            "tf": "CTCF, public reference",
            "reference_type": "chipseq",
            "peak_file": str(peaks.resolve()),
        }
    ]


def test_public_chipseq_validation_does_not_load_full_fasta_sequences(
    tmp_path,
    monkeypatch,
    capsys,
):
    import peak_enrichment

    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGT\nACGT\n", encoding="utf-8")
    peaks = tmp_path / "reference.bed"
    peaks.write_text("chr1\t1\t7\n", encoding="utf-8")
    manifest = tmp_path / "chipseq.csv"
    manifest.write_text(
        "reference_id,tf,peak_file\npublic_ctcf,CTCF,reference.bed\n",
        encoding="utf-8",
    )

    def reject_full_sequence_loading(*args, **kwargs):
        raise AssertionError("preflight validation must use chromosome lengths only")

    monkeypatch.setattr(
        peak_enrichment,
        "load_fasta_sequences",
        reject_full_sequence_loading,
    )

    exit_code = peak_enrichment.main(
        [
            "validate-chipseq-manifest",
            "--manifest",
            str(manifest),
            "--fasta",
            str(fasta),
        ]
    )

    assert exit_code == 0
    assert capsys.readouterr().err == ""


def test_load_reference_manifest_rejects_duplicate_ids_and_invalid_peak_rows(tmp_path):
    from peak_enrichment import load_reference_manifest

    first = tmp_path / "first.bed"
    first.write_text("chr1\t10\t20\n", encoding="utf-8")
    second = tmp_path / "second.bed"
    second.write_text("chrX\t10\t20\n", encoding="utf-8")

    duplicate_manifest = tmp_path / "duplicate.tsv"
    duplicate_manifest.write_text(
        "reference_id\ttf\tpeak_file\n"
        "dup\tCTCF\tfirst.bed\n"
        "dup\tGATA1\tfirst.bed\n",
        encoding="utf-8",
    )
    invalid_peak_manifest = tmp_path / "invalid.tsv"
    invalid_peak_manifest.write_text(
        "reference_id\ttf\tpeak_file\n"
        "bad\tGATA1\tsecond.bed\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate reference_id"):
        load_reference_manifest(duplicate_manifest, {"chr1": 100})

    with pytest.raises(ValueError, match="unknown chromosome"):
        load_reference_manifest(invalid_peak_manifest, {"chr1": 100})


def test_cli_writes_all_four_matrices_plots_and_status_output(tmp_path):
    outdir = run_enrichment_fixture(tmp_path)

    for model in ("random", "length_matched", "gc_matched", "length_gc_matched"):
        assert (outdir / f"matrix_{model}.tsv").exists()
        assert (outdir / f"matrix_{model}.png").exists()
    assert (outdir / "peak_enrichment.tsv").exists()
    assert (outdir / "observed_vs_null.png").exists()
    assert (outdir / "enrichment_status.tsv").exists()

    rows = _read_tsv(outdir / "peak_enrichment.tsv")
    assert len(rows) == 4
    assert list(rows[0]) == [
        "foreground_id",
        "foreground_tf",
        "reference_id",
        "reference_type",
        "reference_tf",
        "background_model",
        "foreground_peak_count",
        "reference_peak_count",
        "observed_overlap_count",
        "null_mean_overlap",
        "null_sd_overlap",
        "enrichment_ratio",
        "empirical_p_value",
        "permutations_requested",
        "permutations_succeeded",
        "seed",
        "status",
    ]
    assert {row["reference_type"] for row in rows} == {"chipseq"}
    assert {row["foreground_peak_count"] for row in rows} == {"2"}
    assert {row["reference_peak_count"] for row in rows} == {"2"}
    assert {row["reference_id"] for row in rows} == {"ref_ctcf_prior"}
    assert {row["background_model"] for row in rows} == {
        "random",
        "length_matched",
        "gc_matched",
        "length_gc_matched",
    }

    matrix_rows = _read_tsv(outdir / "matrix_random.tsv")
    assert matrix_rows[0]["foreground_id"] == "fg_ctcf"
    assert "ref_ctcf_prior" in matrix_rows[0]
    assert "fg_ctcf" not in matrix_rows[0]


def test_observed_vs_null_plot_draws_null_sd_uncertainty():
    from matplotlib.container import ErrorbarContainer

    from peak_enrichment import _load_pyplot, _plot_observed_vs_null

    pyplot = _load_pyplot()
    figure, axis = pyplot.subplots()
    _plot_observed_vs_null(
        axis,
        [
            {
                "foreground_id": "fg",
                "reference_id": "ref",
                "background_model": "random",
                "status": "ok",
                "observed_overlap_count": 4,
                "null_mean_overlap": 2.0,
                "null_sd_overlap": 0.5,
            }
        ],
    )

    assert any(isinstance(container, ErrorbarContainer) for container in axis.containers)
    pyplot.close(figure)


def test_cli_emits_status_only_outputs_when_foreground_manifest_has_no_rows(tmp_path):
    from peak_enrichment import main

    outdir = tmp_path / "enrichment"
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + ("ACGT" * 50) + "\n", encoding="utf-8")

    foreground_manifest = tmp_path / "foregrounds.tsv"
    foreground_manifest.write_text(
        "foreground_id\tforeground_tf\tpeak_file\n",
        encoding="utf-8",
    )

    reference_dir = tmp_path / "references"
    reference_dir.mkdir()
    (reference_dir / "ctcf_prior.bed").write_text(
        "chr1\t12\t18\nchr1\t70\t80\n",
        encoding="utf-8",
    )
    reference_manifest = tmp_path / "references.tsv"
    reference_manifest.write_text(
        "reference_id\ttf\treference_type\tpeak_file\n"
        "ref_ctcf_prior\tCTCF\tchipseq\treferences/ctcf_prior.bed\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "--foreground-manifest",
            str(foreground_manifest),
            "--reference-manifest",
            str(reference_manifest),
            "--fasta",
            str(fasta),
            "--outdir",
            str(outdir),
            "--permutations",
            "8",
            "--seed",
            "1729",
            "--gc-tolerance",
            "0.02",
        ]
    )

    assert exit_code == 0
    assert _read_tsv(outdir / "peak_enrichment.tsv") == []
    status_rows = _read_tsv(outdir / "enrichment_status.tsv")
    assert len(status_rows) == 4
    assert {row["status"] for row in status_rows} == {"no_foreground_peaks"}
    assert {row["background_model"] for row in status_rows} == {
        "random",
        "length_matched",
        "gc_matched",
        "length_gc_matched",
    }
    for model in ("random", "length_matched", "gc_matched", "length_gc_matched"):
        assert (outdir / f"matrix_{model}.tsv").exists()
        assert (outdir / f"matrix_{model}.png").exists()
    assert (outdir / "observed_vs_null.png").exists()
