import csv
import contextlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))
GOLDEN_DEEPTOOLS_PROFILE = (
    ROOT / "tests" / "data" / "qc" / "deeptools_3.5.5_plotprofile.tsv"
)

import qc_dashboard as qc


def library_metrics_tsv(sample_id, **overrides):
    values = {
        "raw_total_reads": 100,
        "mapped_percent": 95,
        "properly_paired_percent": 90,
        "mapq_filtered_reads": 80,
        "mapq_filtered_fragments": 40,
        "mapq_filtered_fraction": 0.8,
        "markdup_examined_reads": 100,
        "duplicate_total": 10,
        "duplicate_percent": 10,
        "mitochondrial_percent": 2,
        "estimated_library_size": 1000,
        "insert_size_total_pairs": 1,
        "insert_size_min": 147,
        "insert_size_q25": 147,
        "insert_size_mean": 147,
        "insert_size_median": 147,
        "insert_size_q75": 147,
        "insert_size_max": 147,
    }
    values.update(overrides)
    columns = ["sample_id", *values]
    return (
        "\t".join(columns) + "\n"
        + "\t".join(
            str(sample_id if column == "sample_id" else values[column])
            for column in columns
        )
        + "\n"
    )


def peak_metrics_tsv(sample_id, **overrides):
    values = {
        "peak_count": 1,
        "total_covered_bases": 321,
        "peak_width_min": 321,
        "peak_width_mean": 321,
        "peak_width_median": 321,
        "peak_width_max": 321,
        "peak_width_q25": 321,
        "peak_width_q75": 321,
        "peak_score_count": 1,
        "peak_score_min": 10,
        "peak_score_q25": 10,
        "peak_score_mean": 10,
        "peak_score_median": 10,
        "peak_score_q75": 10,
        "peak_score_max": 10,
        "signal_value_count": 1,
        "signal_value_min": 5,
        "signal_value_q25": 5,
        "signal_value_mean": 5,
        "signal_value_median": 5,
        "signal_value_q75": 5,
        "signal_value_max": 5,
        "total_fragments": 4,
        "fragments_in_peaks": 1,
        "frip": 0.25,
    }
    values.update(overrides)
    return (
        "metric\tvalue\n"
        f"sample_id\t{sample_id}\n"
        + "".join(f"{metric}\t{value}\n" for metric, value in values.items())
    )


class MetadataAndTableParserTests(unittest.TestCase):
    def make_workspace(self):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        return Path(temporary_directory.name)

    def test_load_metadata_rejects_duplicate_sample_ids(self):
        """A duplicate derived sample would make all later joins ambiguous."""
        workspace = self.make_workspace()
        path = workspace / "metadata.json"
        path.write_text(json.dumps([
            {"sample_id": "S1", "library_id": "L1", "assay_target": "CTCF", "is_control": False},
            {"sample_id": "S1", "library_id": "L1", "assay_target": "CTCF", "is_control": False},
        ]), encoding="utf-8")

        with self.assertRaisesRegex(qc.DashboardInputError, "duplicate sample_id S1"):
            qc.load_metadata(path)

    def test_load_metadata_rejects_blank_required_identity(self):
        """Blank sample identities cannot be used as report-data join keys."""
        workspace = self.make_workspace()
        path = workspace / "metadata.json"
        path.write_text(json.dumps([
            {"sample_id": "", "library_id": "L1", "assay_target": "CTCF", "is_control": False},
        ]), encoding="utf-8")

        with self.assertRaisesRegex(qc.DashboardInputError, "sample_id is required"):
            qc.load_metadata(path)

    def test_load_metadata_rejects_contradictory_control_and_target_links(self):
        """Control/target roles must form a same-input, unambiguous biological link."""
        workspace = self.make_workspace()
        cases = (
            (
                "control-not-igg",
                [{"sample_id": "S_IgG", "library_id": "L1", "input_group": "25K",
                  "assay_target": "CTCF", "is_control": True, "control_id": "", "expected_motif": ""}],
                "control S_IgG assay_target must be IgG",
            ),
            (
                "control-has-target-fields",
                [{"sample_id": "S_IgG", "library_id": "L1", "input_group": "25K",
                  "assay_target": "IgG", "is_control": True, "control_id": "S_IgG", "expected_motif": "CTCF"}],
                "control S_IgG control_id must be blank",
            ),
            (
                "control-has-expected-motif",
                [{"sample_id": "S_IgG", "library_id": "L1", "input_group": "25K",
                  "assay_target": "IgG", "is_control": True, "control_id": "", "expected_motif": "CTCF"}],
                "control S_IgG expected_motif must be blank",
            ),
            (
                "target-missing-control",
                [{"sample_id": "S1", "library_id": "L1", "input_group": "25K",
                  "assay_target": "CTCF", "is_control": False, "control_id": "", "expected_motif": "CTCF"}],
                "target S1 control_id is required",
            ),
            (
                "target-missing-expected-motif",
                [
                    {"sample_id": "S_IgG", "library_id": "L1", "input_group": "25K",
                     "assay_target": "IgG", "is_control": True, "control_id": "", "expected_motif": ""},
                    {"sample_id": "S1", "library_id": "L1", "input_group": "25K",
                     "assay_target": "CTCF", "is_control": False, "control_id": "S_IgG", "expected_motif": ""},
                ],
                "target S1 expected_motif is required",
            ),
            (
                "target-cross-input-control",
                [
                    {"sample_id": "S_IgG", "library_id": "L1", "input_group": "50K",
                     "assay_target": "IgG", "is_control": True, "control_id": "", "expected_motif": ""},
                    {"sample_id": "S1", "library_id": "L2", "input_group": "25K",
                     "assay_target": "CTCF", "is_control": False, "control_id": "S_IgG", "expected_motif": "CTCF"},
                ],
                "target S1 control_id S_IgG must share input_group",
            ),
        )
        for name, metadata, expression in cases:
            with self.subTest(name=name):
                path = workspace / f"{name}.json"
                path.write_text(json.dumps(metadata), encoding="utf-8")
                with self.assertRaisesRegex(qc.DashboardInputError, expression):
                    qc.load_metadata(path)

    def metadata(self):
        return {
            "S1": {
                "sample_id": "S1", "library_id": "L1", "input_group": "25K",
                "assay_target": "CTCF", "is_control": False, "control_id": "S_IgG",
                "expected_motif": "CTCF",
            },
            "S2": {
                "sample_id": "S2", "library_id": "L1", "input_group": "25K",
                "assay_target": "GATA1", "is_control": False, "control_id": "S_IgG",
                "expected_motif": "GATA1",
            },
            "S_IgG": {
                "sample_id": "S_IgG", "library_id": "L1", "input_group": "25K",
                "assay_target": "IgG", "is_control": True, "control_id": None,
                "expected_motif": None,
            },
        }

    def test_read_demultiplex_metrics_preserves_library_and_sample_counts(self):
        """Cached fractions cannot override the physical-library count totals."""
        workspace = self.make_workspace()
        path = workspace / "demultiplex.metrics.json"
        path.write_text(json.dumps({
            "library_id": "L1",
            "total_reads": 100,
            "assigned_reads": 80,
            "ambiguous_reads": 5,
            "unassigned_reads": 15,
            "assigned_fraction": 0.01,
            "assignment_counts": {"S1": 30, "S2": 50, "S_IgG": 0},
        }), encoding="utf-8")

        rows = qc.read_demultiplex_metrics([path], self.metadata())

        self.assertEqual(rows["L1"]["assigned_fraction"], 0.8)
        self.assertEqual(rows["L1"]["assignment_counts"]["S1"], 30)

    def test_read_demultiplex_metrics_preserves_one_argument_embedded_identity_api(self):
        """The declared public parser accepts legacy JSON carrying library_id."""
        workspace = self.make_workspace()
        path = workspace / "demultiplex.metrics.json"
        path.write_text(json.dumps({
            "library_id": "L1", "total_reads": 100, "assigned_reads": 80,
            "ambiguous_reads": 5, "unassigned_reads": 15,
            "assignment_counts": {"S1": 30, "S2": 50},
        }), encoding="utf-8")

        rows = qc.read_demultiplex_metrics([path])

        self.assertEqual(rows["L1"]["assigned_fraction"], 0.8)

    def test_read_demultiplex_metrics_requires_metadata_for_producer_schema_without_identity(self):
        """Producer JSON lacks library_id, so its identity cannot be inferred alone."""
        workspace = self.make_workspace()
        path = workspace / "demultiplex.metrics.json"
        path.write_text(json.dumps({
            "total_reads": 100, "assigned_reads": 80,
            "ambiguous_reads": 5, "unassigned_reads": 15,
            "assignment_counts": {"S1": 30, "S2": 50},
        }), encoding="utf-8")

        with self.assertRaisesRegex(
            qc.DashboardInputError, "validated metadata is required to infer library_id"
        ):
            qc.read_demultiplex_metrics([path])

    def test_read_demultiplex_metrics_rejects_duplicate_library_ids(self):
        """A library must contribute exactly one unambiguous count record."""
        workspace = self.make_workspace()
        paths = []
        for index in range(2):
            path = workspace / f"demultiplex-{index}.json"
            path.write_text(json.dumps({
                "library_id": "L1", "total_reads": 1, "assigned_reads": 1,
                "ambiguous_reads": 0, "unassigned_reads": 0,
                "assignment_counts": {"S1": 1, "S2": 0, "S_IgG": 0},
            }), encoding="utf-8")
            paths.append(path)

        with self.assertRaisesRegex(qc.DashboardInputError, "duplicate library_id L1"):
            qc.read_demultiplex_metrics(paths, self.metadata())

    def test_read_demultiplex_metrics_rejects_negative_or_nonfinite_counts(self):
        """Counts must remain finite non-negative quantities in fractions and joins."""
        workspace = self.make_workspace()
        for name, total_reads, expression in (
            ("negative", -1, "finite and >= 0"),
            ("nonfinite", "nan", "finite and >= 0"),
        ):
            with self.subTest(name=name):
                path = workspace / f"{name}.json"
                path.write_text(json.dumps({
                    "library_id": "L1", "total_reads": total_reads, "assigned_reads": 0,
                    "ambiguous_reads": 0, "unassigned_reads": 0, "assignment_counts": {},
                }), encoding="utf-8")
                with self.assertRaisesRegex(qc.DashboardInputError, expression):
                    qc.read_demultiplex_metrics([path], self.metadata())

    def test_read_demultiplex_metrics_recovers_library_id_from_real_producer_schema(self):
        """The producer JSON has only count fields and per-sample assignments."""
        workspace = self.make_workspace()
        path = workspace / "demultiplex.metrics.json"
        path.write_text(json.dumps({
            "total_reads": 100,
            "assigned_reads": 80,
            "ambiguous_reads": 5,
            "unassigned_reads": 15,
            "assignment_counts": {"S1": 30, "S2": 50, "S_IgG": 0},
            "observed_i2_counts": {"TATAGCCT": 30, "ATAGAGGC": 50},
        }), encoding="utf-8")

        rows = qc.read_demultiplex_metrics([path], self.metadata())

        self.assertEqual(
            rows["L1"]["assignment_counts"],
            {"S1": 30, "S2": 50, "S_IgG": 0},
        )
        self.assertEqual(rows["L1"]["assigned_fraction"], 0.8)

    def test_read_demultiplex_metrics_requires_every_metadata_sample_allocation(self):
        """A dropped zero-count sample key must not be reported as computed."""
        workspace = self.make_workspace()
        path = workspace / "demultiplex.metrics.json"
        path.write_text(json.dumps({
            "total_reads": 100,
            "assigned_reads": 80,
            "ambiguous_reads": 5,
            "unassigned_reads": 15,
            "assignment_counts": {"S1": 30, "S2": 50},
        }), encoding="utf-8")

        with self.assertRaisesRegex(
            qc.DashboardInputError,
            "assignment_counts for L1 must exactly match metadata samples",
        ):
            qc.read_demultiplex_metrics([path], self.metadata())

    def test_read_demultiplex_metrics_rejects_ambiguous_unknown_or_inconsistent_counts(self):
        """Assignment identities and count totals must be safe before joining samples."""
        workspace = self.make_workspace()
        metadata = self.metadata()
        metadata["S2"] = {**metadata["S2"], "library_id": "L2"}
        cases = (
            (
                "ambiguous-library",
                metadata,
                {"S1": 30, "S2": 50},
                80, 5, 15,
                "assignment_counts map to multiple library_id values",
            ),
            (
                "unknown-sample",
                self.metadata(),
                {"S1": 30, "UNKNOWN": 50},
                80, 5, 15,
                "assignment_counts reference unknown sample_id UNKNOWN",
            ),
            (
                "bad-assigned-sum",
                self.metadata(),
                {"S1": 30, "S2": 40},
                80, 5, 15,
                "assignment_counts must sum to assigned_reads",
            ),
            (
                "bad-total-sum",
                self.metadata(),
                {"S1": 30, "S2": 50},
                80, 5, 14,
                "assigned_reads, ambiguous_reads, and unassigned_reads must sum to total_reads",
            ),
        )
        for name, case_metadata, assignments, assigned, ambiguous, unassigned, expression in cases:
            with self.subTest(name=name):
                path = workspace / f"{name}.json"
                path.write_text(json.dumps({
                    "total_reads": 100, "assigned_reads": assigned,
                    "ambiguous_reads": ambiguous, "unassigned_reads": unassigned,
                    "assignment_counts": assignments,
                }), encoding="utf-8")
                with self.assertRaisesRegex(qc.DashboardInputError, expression):
                    qc.read_demultiplex_metrics([path], case_metadata)

    def test_read_library_metrics_rejects_duplicate_sample_ids(self):
        """A sample may have only one library-QC summary row."""
        workspace = self.make_workspace()
        paths = []
        for index in range(2):
            path = workspace / f"library-{index}.tsv"
            path.write_text(library_metrics_tsv("S1"), encoding="utf-8")
            paths.append(path)

        with self.assertRaisesRegex(qc.DashboardInputError, "duplicate sample_id S1"):
            qc.read_library_metrics(paths)

    def test_read_peak_metrics_parses_frip_and_rejects_duplicate_metric_names(self):
        """Metric/value peak tables require one value per metric name."""
        workspace = self.make_workspace()
        valid = workspace / "S1.peak_qc.tsv"
        valid.write_text(
            peak_metrics_tsv("S1", peak_count=12, frip=0.42),
            encoding="utf-8",
        )
        self.assertEqual(qc.read_peak_metrics([valid])["S1"]["frip"], 0.42)
        self.assertEqual(qc.read_peak_metrics([valid])["S1"]["width_min"], 321)
        self.assertNotIn("peak_width_min", qc.read_peak_metrics([valid])["S1"])

        duplicate = workspace / "duplicate.peak_qc.tsv"
        duplicate.write_text(
            peak_metrics_tsv("S2") + "frip\t0.43\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(qc.DashboardInputError, "duplicate peak metric frip"):
            qc.read_peak_metrics([duplicate])

    def test_library_and_peak_parsers_reject_incomplete_producer_schemas(self):
        """A present but structurally empty artifact must not become computed QC."""
        workspace = self.make_workspace()
        library = workspace / "S1.library_qc.tsv"
        library.write_text("sample_id\nS1\n", encoding="utf-8")
        peak = workspace / "S1.peak_qc.tsv"
        peak.write_text(
            "metric\tvalue\nsample_id\tS1\npeak_count\t1\nfrip\t0.25\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(
            qc.DashboardInputError,
            "library metrics columns must exactly match the producer schema",
        ):
            qc.read_library_metrics([library])
        with self.assertRaisesRegex(
            qc.DashboardInputError,
            "peak metrics are missing required producer fields",
        ):
            qc.read_peak_metrics([peak])

    def test_distribution_parsers_consume_real_producer_tables_and_preserve_empty_samples(self):
        """Dropping either staged distribution would remove approved report data."""
        workspace = self.make_workspace()
        insert = workspace / "S1.insert_size_distribution.tsv"
        insert.write_text(
            "sample_id\tinsert_size\tpair_count\n"
            "S1\t100\t2\n"
            "S1\t200\t1\n",
            encoding="utf-8",
        )
        empty_insert = workspace / "S2.insert_size_distribution.tsv"
        empty_insert.write_text(
            "sample_id\tinsert_size\tpair_count\n",
            encoding="utf-8",
        )
        widths = workspace / "S1.peak_qc.width_histogram.tsv"
        widths.write_text(
            "width\tpeak_count\n"
            "150\t3\n"
            "300\t1\n",
            encoding="utf-8",
        )
        empty_widths = workspace / "S2.peak_qc.width_histogram.tsv"
        empty_widths.write_text("width\tpeak_count\n", encoding="utf-8")

        self.assertEqual(
            qc.read_insert_size_distributions([insert, empty_insert]),
            {
                "S1": [
                    {"insert_size": 100, "pair_count": 2},
                    {"insert_size": 200, "pair_count": 1},
                ],
                "S2": [],
            },
        )
        self.assertEqual(
            qc.read_peak_width_distributions([widths, empty_widths]),
            {
                "S1": [
                    {"width": 150, "peak_count": 3},
                    {"width": 300, "peak_count": 1},
                ],
                "S2": [],
            },
        )

    def test_fragments_per_peak_parser_compacts_real_producer_rows(self):
        """Raw per-peak identities must collapse to a small count histogram."""
        workspace = self.make_workspace()
        path = workspace / "S1.peak_qc.fragments_per_peak.tsv"
        path.write_text(
            "chrom\tstart\tend\tpeak_name\twidth\tscore\tsignal_value\tfragment_count\n"
            "chr1\t0\t250\tp1\t250\t10\t3.5\t0\n"
            "chr1\t500\t1000\tp2\t500\t12\t5.0\t4\n"
            "chr2\t0\t750\tp3\t750\t8\t2.0\t4\n",
            encoding="utf-8",
        )
        empty_path = workspace / "S2.peak_qc.fragments_per_peak.tsv"
        empty_path.write_text(
            "chrom\tstart\tend\tpeak_name\twidth\tscore\t"
            "signal_value\tfragment_count\n",
            encoding="utf-8",
        )

        self.assertEqual(
            qc.read_fragments_per_peak_distributions([path, empty_path]),
            {
                "S1": [
                    {"fragment_count": 0, "peak_count": 1},
                    {"fragment_count": 4, "peak_count": 2},
                ],
                "S2": [],
            },
        )

    def test_fragments_per_peak_parser_rejects_ambiguous_or_invalid_inputs(self):
        """Malformed producer data must not be silently joined to a target."""
        workspace = self.make_workspace()
        header = (
            "chrom\tstart\tend\tpeak_name\twidth\tscore\t"
            "signal_value\tfragment_count\n"
        )
        cases = (
            (
                "wrong-header.peak_qc.fragments_per_peak.tsv",
                "chrom\tstart\tend\tpeak_name\tfragment_count\n",
                "columns must exactly match the producer schema",
            ),
            (
                "negative.peak_qc.fragments_per_peak.tsv",
                header + "chr1\t0\t1\tp1\t1\t1\t1\t-1\n",
                "finite and >= 0",
            ),
            (
                "fractional.peak_qc.fragments_per_peak.tsv",
                header + "chr1\t0\t1\tp1\t1\t1\t1\t1.5\n",
                "must be an integer",
            ),
            (
                "wrong.tsv",
                header,
                "unexpected fragments-per-peak distribution filename",
            ),
        )
        for filename, contents, expression in cases:
            with self.subTest(filename=filename):
                path = workspace / filename
                path.write_text(contents, encoding="utf-8")
                with self.assertRaisesRegex(qc.DashboardInputError, expression):
                    qc.read_fragments_per_peak_distributions([path])

        duplicate_a = workspace / "first" / "S1.peak_qc.fragments_per_peak.tsv"
        duplicate_b = workspace / "second" / "S1.peak_qc.fragments_per_peak.tsv"
        duplicate_a.parent.mkdir()
        duplicate_b.parent.mkdir()
        duplicate_a.write_text(header, encoding="utf-8")
        duplicate_b.write_text(header, encoding="utf-8")
        with self.assertRaisesRegex(
            qc.DashboardInputError,
            "duplicate fragments-per-peak distribution sample_id S1",
        ):
            qc.read_fragments_per_peak_distributions([duplicate_a, duplicate_b])

    def test_distribution_parsers_reject_filename_and_row_identity_mismatch(self):
        """A row must never be attached to the sample encoded by another file."""
        workspace = self.make_workspace()
        path = workspace / "S1.insert_size_distribution.tsv"
        path.write_text(
            "sample_id\tinsert_size\tpair_count\nS2\t100\t2\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(
            qc.DashboardInputError,
            "distribution sample_id S2 does not match filename sample_id S1",
        ):
            qc.read_insert_size_distributions([path])


class TssAndAmeTests(unittest.TestCase):
    def make_workspace(self):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        return Path(temporary_directory.name)

    def test_tss_enrichment_uses_center_and_terminal_100bp(self):
        """The pinned deepTools table drives the exact center/flank calculation."""
        profile = qc.read_tss_profile(GOLDEN_DEEPTOOLS_PROFILE)

        self.assertEqual(profile[0], (-3000, 2.0))
        self.assertEqual(profile[300], (0, 12.0))
        self.assertEqual(qc.calculate_tss_enrichment(profile), 6.0)

    def test_read_tss_profile_accepts_deeptools_decimal_bin_axis(self):
        """deepTools 3.5.5 serializes its NumPy-generated bin axis as floats."""
        workspace = self.make_workspace()
        path = workspace / "decimal-bins.tsv"
        lines = GOLDEN_DEEPTOOLS_PROFILE.read_text(encoding="utf-8").splitlines()
        bins = lines[1].split("\t")
        bins[2:] = [f"{value}.0" for value in bins[2:]]
        path.write_text(
            "\n".join((lines[0], "\t".join(bins), lines[2])) + "\n",
            encoding="utf-8",
        )

        profile = qc.read_tss_profile(path)

        self.assertEqual(profile[0], (-3000, 2.0))
        self.assertEqual(profile[300], (0, 12.0))
        self.assertEqual(profile[-1], (2990, 2.0))

    def test_tss_zero_flank_returns_missing(self):
        """A zero baseline has no defined center-to-flank enrichment ratio."""
        profile = [(index * 10 - 3000, 0.0) for index in range(600)]

        self.assertIsNone(qc.calculate_tss_enrichment(profile))

    def test_read_tss_profile_rejects_wrong_bin_count_and_noncontiguous_axis(self):
        """The real deepTools header and one-based axis are part of the contract."""
        workspace = self.make_workspace()
        wrong_count = workspace / "wrong-count.tsv"
        lines = GOLDEN_DEEPTOOLS_PROFILE.read_text(encoding="utf-8").splitlines()
        wrong_count.write_text(
            "\n".join("\t".join(line.split("\t")[:-1]) for line in lines) + "\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(qc.DashboardInputError, "602 tab-separated fields"):
            qc.read_tss_profile(wrong_count)

        noncontiguous = workspace / "noncontiguous.tsv"
        bins = lines[1].split("\t")
        bins[302] = "999"
        noncontiguous.write_text(
            "\n".join((lines[0], "\t".join(bins), lines[2])) + "\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(qc.DashboardInputError, "bins must be exactly 1 through 600"):
            qc.read_tss_profile(noncontiguous)

    def test_read_tss_profile_rejects_missing_or_multiple_profiles(self):
        """The dashboard needs exactly one aggregate profile row per input file."""
        workspace = self.make_workspace()
        missing = workspace / "missing.tsv"
        missing.write_text("# no aggregate profile\n\n", encoding="utf-8")
        with self.assertRaisesRegex(qc.DashboardInputError, "exactly three non-comment rows"):
            qc.read_tss_profile(missing)

        multiple = workspace / "multiple.tsv"
        lines = GOLDEN_DEEPTOOLS_PROFILE.read_text(encoding="utf-8").splitlines()
        multiple.write_text("\n".join((*lines, lines[2])) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(qc.DashboardInputError, "exactly three non-comment rows"):
            qc.read_tss_profile(multiple)

    def test_read_tss_profile_rejects_raw_multiregion_matrix(self):
        """Only Task 3's aggregate plotProfile table is a supported dashboard input."""
        workspace = self.make_workspace()
        path = workspace / "S1.tss_matrix.tsv"
        values = "\t".join(["1"] * 600)
        path.write_text(
            "# raw computeMatrix rows are per-region rather than aggregate\n"
            f"chr1\t0\t1\tregion_1\t0\t+\t{values}\n"
            f"chr1\t2\t3\tregion_2\t0\t+\t{values}\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(qc.DashboardInputError, "exactly three non-comment rows"):
            qc.read_tss_profile(path)

    def test_nonfinite_and_zero_tss_flanks_become_missing_with_diagnostics(self):
        """Unusable flank signal is reportable missing data, not a fatal join error."""
        workspace = self.make_workspace()
        lines = GOLDEN_DEEPTOOLS_PROFILE.read_text(encoding="utf-8").splitlines()
        for sample_id, replacement, expected_status in (
            ("ZERO", "0", "zero_flank"),
            ("NONFINITE", "nan", "non_finite_flank"),
        ):
            data = lines[2].split("\t")
            data[0] = "coverage"
            data[2:12] = [replacement] * 10
            data[-10:] = [replacement] * 10
            (workspace / f"{sample_id}.tss_profile.tsv").write_text(
                "\n".join((lines[0], lines[1], "\t".join(data))) + "\n",
                encoding="utf-8",
            )
            (workspace / f"{sample_id}.tss_status.tsv").write_text(
                "sample_id\tannotation_mode\tstatus\n"
                f"{sample_id}\tbed\tcomputed\n",
                encoding="utf-8",
            )

        metrics = qc._read_tss_metrics(workspace)

        self.assertIsNone(metrics["ZERO"]["enrichment"])
        self.assertEqual(metrics["ZERO"]["score_status"], "zero_flank")
        self.assertIsNone(metrics["NONFINITE"]["enrichment"])
        self.assertEqual(metrics["NONFINITE"]["score_status"], "non_finite_flank")

    def test_tss_status_profile_identity_mismatches_are_explicit_failed_records(self):
        """Partial TSS artifacts must not be silently promoted to computed results."""
        workspace = self.make_workspace()
        (workspace / "STATUS_ONLY.tss_status.tsv").write_text(
            "sample_id\tannotation_mode\tstatus\n"
            "STATUS_ONLY\tbed\tcomputed\n",
            encoding="utf-8",
        )
        shutil.copyfile(
            GOLDEN_DEEPTOOLS_PROFILE,
            workspace / "PROFILE_ONLY.tss_profile.tsv",
        )

        metrics = qc._read_tss_metrics(workspace)

        self.assertEqual(metrics["STATUS_ONLY"]["status"], "failed")
        self.assertEqual(metrics["STATUS_ONLY"]["reason"], "profile_missing")
        self.assertEqual(metrics["PROFILE_ONLY"]["status"], "failed")
        self.assertEqual(metrics["PROFILE_ONLY"]["reason"], "status_missing")

    def test_fake_plotprofile_output_is_byte_identical_to_pinned_golden_fixture(self):
        """The offline tool must not mask production-format parser defects."""
        workspace = self.make_workspace()
        profile = workspace / "profile.png"
        table = workspace / "profile.tsv"

        completed = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tests" / "data" / "e2e" / "fakebin" / "fake_bio_tool.py"),
                "--fake-tool", "plotProfile",
                "--matrixFile", str(workspace / "matrix.gz"),
                "--outFileName", str(profile),
                "--outFileNameData", str(table),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(table.read_bytes(), GOLDEN_DEEPTOOLS_PROFILE.read_bytes())

    def test_read_top_ame_ignores_no_peaks_and_orders_top_ten(self):
        """AME rankings must be deterministic and omit the pipeline's no-peaks sentinel."""
        workspace = self.make_workspace()
        path = workspace / "ame.tsv"
        rows = [
            "rank\tmotif_ID\tmotif_Alt_ID\tadj_p-value\tp-value\tscore\tpos",
            "1\t__NO_PEAKS__\t__NO_PEAKS__\t0\t0\t0\t0",
            "2\tB\tb\t0.001\t0.001\t2\t4",
            "2\tA\ta\t0.001\t0.001\t2\t4",
        ]
        rows.extend(
            f"{rank}\tM{rank}\tm{rank}\t0.01\t0.01\t2\t4"
            for rank in range(1, 11)
        )
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")

        top = qc.read_top_ame(path)

        self.assertEqual(len(top), 10)
        self.assertEqual([row["motif_id"] for row in top[:2]], ["A", "B"])
        self.assertNotIn("__NO_PEAKS__", [row["motif_id"] for row in top])
        self.assertEqual(
            set(top[0]),
            {"motif_id", "motif_alt_id", "adjusted_p_value", "p_value", "effect", "positive_sequences", "rank"},
        )

    def test_read_all_ame_retains_rows_beyond_top_ten(self):
        """The heatmap needs the complete AME result while legacy exports stay bounded."""
        workspace = self.make_workspace()
        path = workspace / "ame.tsv"
        rows = [
            "rank\tmotif_ID\tmotif_Alt_ID\tadj_p-value\tp-value\tscore\tpos",
            "1\t__NO_PEAKS__\t__NO_PEAKS__\t0\t0\t0\t0",
        ]
        rows.extend(
            f"{rank}\tM{rank:02d}\tmotif-{rank:02d}\t{rank / 1000}\t"
            f"{rank / 100}\t2\t4"
            for rank in range(1, 13)
        )
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")

        all_records = qc.read_all_ame(path)
        top_records = qc.read_top_ame(path)

        self.assertEqual(len(all_records), 12)
        self.assertEqual(len(top_records), 10)
        self.assertEqual(
            [row["motif_id"] for row in top_records],
            [row["motif_id"] for row in all_records[:10]],
        )
        self.assertNotIn("__NO_PEAKS__", [row["motif_id"] for row in all_records])

    def test_cognate_matching_uses_complete_case_insensitive_tokens(self):
        """Punctuation-delimited TF names match, but longer symbols do not."""
        self.assertTrue(qc.is_cognate_motif("GATA1", "MA0140.2", "GATA1::TAL1"))
        self.assertTrue(qc.is_cognate_motif("GATA1", "TAL1::gata1", "complex"))
        self.assertFalse(qc.is_cognate_motif("GATA1", "MA9999", "GATA10"))

    def test_ame_directory_load_keeps_complete_internal_and_bounded_public_records(self):
        """A single AME read supplies both heatmap data and the legacy top-ten export."""
        workspace = self.make_workspace()
        path = workspace / "S1" / "ame" / "ame.tsv"
        path.parent.mkdir(parents=True)
        rows = [
            "rank\tmotif_ID\tmotif_Alt_ID\tadj_p-value\tp-value\tscore\tpos",
        ]
        rows.extend(
            f"{rank}\tM{rank:02d}\tmotif-{rank:02d}\t{rank / 1000}\t"
            f"{rank / 100}\t2\t4"
            for rank in range(1, 13)
        )
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")

        all_motifs, top_motifs = qc._read_ame_motifs_from_directory(workspace)
        data = qc.build_report_data(
            {
                "S1": {
                    "sample_id": "S1", "library_id": "L1",
                    "assay_target": "GATA1", "is_control": False,
                    "input_group": None, "control_id": None,
                    "expected_motif": "GATA1",
                },
            },
            {}, {}, {}, {}, {}, top_motifs,
            annotation_status="skipped_no_annotation",
            ame_motifs=all_motifs,
            motif_analysis_status="computed",
        )

        self.assertEqual(len(data["samples_by_id"]["S1"]["ame_motifs"]), 12)
        self.assertEqual(len(data["samples_by_id"]["S1"]["top_motifs"]), 10)
        public_sample = qc._json_payload(data)["samples"][0]
        self.assertNotIn("ame_motifs", public_sample)
        self.assertEqual(len(public_sample["top_motifs"]), 10)

    def test_build_motif_heatmap_selects_targets_cognates_and_top_noncognates(self):
        """The cohort matrix excludes controls and bounds only noncognate rows."""
        noncognate = [
            {
                "motif_id": f"MA_OTHER_{index:02d}",
                "motif_alt_id": f"OTHER{index:02d}",
                "adjusted_p_value": 0.01 + index / 1000 if index < 12 else None,
            }
            for index in range(17)
        ]
        gata_records = [
            {
                "motif_id": "MA_GATA1", "motif_alt_id": "GATA1::TAL1",
                "adjusted_p_value": 0.001,
            },
            {
                "motif_id": "MA_GATA10", "motif_alt_id": "GATA10",
                "adjusted_p_value": 0.15,
            },
            {
                "motif_id": "MA_ZERO", "motif_alt_id": "ZERO",
                "adjusted_p_value": 0,
            },
            {
                "motif_id": "MA_NS", "motif_alt_id": "NS",
                "adjusted_p_value": 0.2,
            },
            *noncognate,
        ]
        samples = [
            {
                "sample_id": "GATA_sample", "expected_motif": "GATA1",
                "is_control": False, "ame_motifs": gata_records,
            },
            {
                "sample_id": "CTCF_sample", "expected_motif": "CTCF",
                "is_control": False,
                "ame_motifs": [{
                    "motif_id": "MA_CTCF", "motif_alt_id": "CTCF",
                    "adjusted_p_value": None,
                }],
            },
            {
                "sample_id": "IgG_sample", "expected_motif": None,
                "is_control": True, "ame_motifs": gata_records,
            },
        ]

        matrix = qc.build_motif_heatmap(samples)

        self.assertNotIn("IgG_sample", matrix["sample_ids"])
        self.assertEqual(matrix["sample_ids"], ["CTCF_sample", "GATA_sample"])
        self.assertIn(("MA_GATA1", "GATA1::TAL1"), matrix["motif_keys"])
        self.assertNotIn(
            ("MA_GATA10", "GATA10"), matrix["forced_cognate_keys"]
        )
        self.assertEqual(len(matrix["noncognate_keys"]), 15)
        self.assertEqual(
            matrix["cells"][("GATA_sample", ("MA_ZERO", "ZERO"))]["score"],
            60.0,
        )
        self.assertEqual(
            matrix["cells"][("GATA_sample", ("MA_ZERO", "ZERO"))]["label"],
            ">60",
        )
        self.assertEqual(
            matrix["cells"][("GATA_sample", ("MA_NS", "NS"))]["label"],
            "ns",
        )
        self.assertEqual(
            matrix["cells"][
                ("GATA_sample", ("MA_GATA1", "GATA1::TAL1"))
            ]["score"],
            3.0,
        )
        self.assertEqual(
            matrix["cells"][("CTCF_sample", ("MA_CTCF", "CTCF"))]["label"],
            "ns",
        )
        self.assertTrue(
            matrix["cells"][
                ("GATA_sample", ("MA_GATA1", "GATA1::TAL1"))
            ]["outlined"]
        )

    def test_build_motif_heatmap_forces_cognate_seen_in_another_sample(self):
        """Every cohort TF forces matching motifs even when its own AME row is absent."""
        samples = [
            {
                "sample_id": "GATA_sample", "expected_motif": "GATA1",
                "is_control": False, "ame_motifs": [],
            },
            {
                "sample_id": "CTCF_sample", "expected_motif": "CTCF",
                "is_control": False,
                "ame_motifs": [{
                    "motif_id": "MA_GATA", "motif_alt_id": "GATA1::TAL1",
                    "adjusted_p_value": 0.001,
                }],
            },
        ]

        matrix = qc.build_motif_heatmap(samples, noncognate_limit=0)
        motif_key = ("MA_GATA", "GATA1::TAL1")

        self.assertEqual(matrix["forced_cognate_keys"], [motif_key])
        self.assertEqual(matrix["motif_keys"], [motif_key])
        self.assertEqual(
            matrix["cells"][("GATA_sample", motif_key)]["label"], "ns"
        )
        self.assertTrue(
            matrix["cells"][("GATA_sample", motif_key)]["outlined"]
        )
        self.assertFalse(
            matrix["cells"][("CTCF_sample", motif_key)]["outlined"]
        )

    def test_build_motif_heatmap_collapses_case_variants_with_stable_display(self):
        """Case variants share one identity while the best value supplies the cell."""
        samples = [{
            "sample_id": "GATA_sample", "expected_motif": "GATA1",
            "is_control": False,
            "ame_motifs": [
                {
                    "motif_id": "ma_gata", "motif_alt_id": "gata1::tal1",
                    "adjusted_p_value": 0.001,
                },
                {
                    "motif_id": "MA_GATA", "motif_alt_id": "GATA1::TAL1",
                    "adjusted_p_value": 0.01,
                },
            ],
        }]

        matrix = qc.build_motif_heatmap(samples, noncognate_limit=0)
        motif_key = ("MA_GATA", "GATA1::TAL1")

        self.assertEqual(matrix["motif_keys"], [motif_key])
        self.assertEqual(matrix["forced_cognate_keys"], [motif_key])
        self.assertEqual(
            matrix["cells"][("GATA_sample", motif_key)]["score"], 3.0
        )


class ReportDataTests(unittest.TestCase):
    def test_build_report_data_joins_target_and_preserves_control_optional_gaps(self):
        """IgG controls retain technical data without being treated as missing biology."""
        metadata = {
            "S_IgG": {
                "sample_id": "S_IgG", "library_id": "L1", "assay_target": "IgG",
                "is_control": True, "input_group": "25K", "control_id": None,
                "expected_motif": None,
            },
            "S_CTCF": {
                "sample_id": "S_CTCF", "library_id": "L1", "assay_target": "CTCF",
                "is_control": False, "input_group": "25K", "control_id": "S_IgG",
                "expected_motif": "CTCF",
            },
        }
        demultiplex = {
            "L1": {
                "total_reads": 100, "assigned_reads": 80, "ambiguous_reads": 5,
                "unassigned_reads": 15, "assigned_fraction": 0.8,
                "ambiguous_fraction": 0.05, "unassigned_fraction": 0.15,
                "assignment_counts": {"S_IgG": 30, "S_CTCF": 50},
            },
        }
        libraries = {"S_CTCF": {"sample_id": "S_CTCF", "mapped_percent": 95.0}}
        peaks = {"S_CTCF": {"sample_id": "S_CTCF", "frip": 0.42}}
        tss = {"S_CTCF": {"status": "computed", "enrichment": 4.2}}
        motifs = {"S_CTCF": {"status": "pass", "best_motif_id": "MA0139.1"}}
        top_motifs = {"S_CTCF": [{"motif_id": "MA0139.1"}]}

        data = qc.build_report_data(
            metadata, demultiplex, libraries, peaks, tss, motifs, top_motifs,
            annotation_status="computed_bed",
        )

        self.assertEqual(data["samples"][0]["sample_id"], "S_CTCF")
        self.assertEqual(
            data["samples_by_id"]["S_CTCF"]["demultiplex"]["sample_assigned_reads"], 50
        )
        self.assertEqual(data["samples_by_id"]["S_CTCF"]["peak"]["frip"], 0.42)
        self.assertIsNone(data["samples_by_id"]["S_IgG"]["peak"]["frip"])
        self.assertEqual(data["samples_by_id"]["S_IgG"]["motif"]["status"], "not_applicable_control")
        self.assertTrue(any("S_IgG" in warning["message"] for warning in data["warnings"]))
        self.assertFalse(any(
            "S_IgG: missing peak metrics" in warning["message"]
            or "S_IgG: missing motif metrics" in warning["message"]
            for warning in data["warnings"]
        ))

    def test_build_report_data_rejects_optional_metrics_for_unknown_samples(self):
        """Optional tables may be absent, but they may not introduce a new identity."""
        metadata = {
            "S1": {
                "sample_id": "S1", "library_id": "L1", "assay_target": "CTCF",
                "is_control": False, "input_group": None, "control_id": None,
                "expected_motif": "CTCF",
            },
        }

        with self.assertRaisesRegex(qc.DashboardInputError, "unknown sample_id S2"):
            qc.build_report_data(
                metadata, {}, {"S2": {"sample_id": "S2"}}, {}, {}, {}, {},
                annotation_status="skipped_no_annotation",
            )
        with self.assertRaisesRegex(qc.DashboardInputError, "unknown sample_id S2"):
            qc.build_report_data(
                metadata, {}, {}, {}, {}, {}, {},
                annotation_status="skipped_no_annotation",
                fragments_per_peak={"S2": []},
            )

    def test_build_report_data_models_intentional_skips_without_generic_missing_status(self):
        """Disabled optional analyses retain a skipped reason instead of looking lost."""
        metadata = {
            "S1": {
                "sample_id": "S1", "library_id": "L1", "assay_target": "CTCF",
                "is_control": False, "input_group": "25K", "control_id": "I1",
                "expected_motif": "CTCF",
            },
            "I1": {
                "sample_id": "I1", "library_id": "L1", "assay_target": "IgG",
                "is_control": True, "input_group": "25K", "control_id": None,
                "expected_motif": None,
            },
        }

        data = qc.build_report_data(
            metadata, {}, {}, {}, {}, {}, {},
            annotation_status="skipped_no_annotation",
            insert_sizes={},
            peak_widths={},
            motif_analysis_status="skipped_no_database",
        )

        target = data["samples_by_id"]["S1"]
        self.assertEqual(
            target["availability"]["tss"],
            {"status": "skipped", "reason": "annotation_not_provided"},
        )
        self.assertEqual(target["tss"]["status"], "skipped_no_annotation")
        self.assertEqual(
            target["availability"]["motif"],
            {"status": "skipped", "reason": "motif_database_not_provided"},
        )
        self.assertFalse(any(
            "missing TSS metrics" in warning["message"]
            or "missing motif metrics" in warning["message"]
            for warning in data["warnings"]
        ))
        self.assertTrue(any(
            warning["family"] == "tss" and warning["status"] == "skipped"
            for warning in data["warnings"]
        ))

    def test_build_report_data_warns_for_partial_empty_failed_and_unusable_optional_results(self):
        """Every unusable optional-family state remains distinguishable in the model."""
        metadata = {
            "A": {
                "sample_id": "A", "library_id": "L1", "assay_target": "CTCF",
                "is_control": False, "input_group": "25K", "control_id": "I",
                "expected_motif": "CTCF",
            },
            "B": {
                "sample_id": "B", "library_id": "L2", "assay_target": "GATA1",
                "is_control": False, "input_group": "25K", "control_id": "I",
                "expected_motif": "GATA1",
            },
            "I": {
                "sample_id": "I", "library_id": "L1", "assay_target": "IgG",
                "is_control": True, "input_group": "25K", "control_id": None,
                "expected_motif": None,
            },
        }
        data = qc.build_report_data(
            metadata,
            {},
            {},
            {"A": {"sample_id": "A", "peak_count": 0, "frip": None}},
            {
                "A": {
                    "sample_id": "A", "status": "computed",
                    "profile": [
                        (index * 10 - 3000, 0.0)
                        for index in range(600)
                    ],
                    "enrichment": None, "score_status": "zero_flank",
                },
                "B": {"sample_id": "B", "status": "failed", "profile": None},
            },
            {
                "A": {"sample_id": "A", "status": "no_peaks", "ame_status": "no_peaks"},
                "B": {"sample_id": "B", "status": "failed", "ame_status": "failed"},
            },
            {"A": [], "B": []},
            annotation_status="computed_bed",
            insert_sizes={"A": [], "B": []},
            peak_widths={"A": []},
            motif_analysis_status="computed",
        )

        a = data["samples_by_id"]["A"]
        b = data["samples_by_id"]["B"]
        self.assertEqual(a["availability"]["insert_size"]["status"], "empty")
        self.assertEqual(a["availability"]["peak_width"]["status"], "empty")
        self.assertEqual(a["availability"]["tss"]["status"], "computed")
        self.assertEqual(a["tss"]["score_status"], "zero_flank")
        self.assertEqual(a["availability"]["motif"]["status"], "empty")
        self.assertEqual(a["availability"]["ame"]["status"], "empty")
        self.assertEqual(b["availability"]["peak_width"]["status"], "missing")
        self.assertEqual(b["availability"]["tss"]["status"], "failed")
        self.assertEqual(b["availability"]["motif"]["status"], "failed")
        self.assertEqual(b["availability"]["ame"]["status"], "failed")
        observed = {
            (warning["sample_id"], warning["family"], warning["status"])
            for warning in data["warnings"]
        }
        self.assertTrue({
            ("A", "insert_size", "empty"),
            ("A", "peak_width", "empty"),
            ("A", "tss", "computed"),
            ("A", "motif", "empty"),
            ("A", "ame", "empty"),
            ("B", "peak_width", "missing"),
            ("B", "tss", "failed"),
            ("B", "motif", "failed"),
            ("B", "ame", "failed"),
        }.issubset(observed))


class DashboardOutputTests(unittest.TestCase):
    """The outputs remain reproducible and safe to open without a network."""

    def make_workspace(self):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        return Path(temporary_directory.name)

    def report_data(self):
        """Hand-authored joined data proves serialization, not parser behavior."""
        metadata = {
            "Z_TARGET": {
                "sample_id": "Z_TARGET", "library_id": "L1", "input_group": "25K",
                "assay_target": "CTCF", "is_control": False, "control_id": "A_IGG",
                "expected_motif": "CTCF & <target>",
            },
            "A_IGG": {
                "sample_id": "A_IGG", "library_id": "L1", "input_group": "25K",
                "assay_target": "IgG", "is_control": True, "control_id": None,
                "expected_motif": None,
            },
        }
        data = qc.build_report_data(
            metadata,
            {
                "L1": {
                    "total_reads": 100, "assigned_reads": 80, "ambiguous_reads": 5,
                    "unassigned_reads": 15, "assigned_fraction": 0.8,
                    "ambiguous_fraction": 0.05, "unassigned_fraction": 0.15,
                    "assignment_counts": {"A_IGG": 30, "Z_TARGET": 50},
                },
            },
            {"Z_TARGET": {"sample_id": "Z_TARGET", "mapped_percent": 95.0}},
            {},
            {
                "A_IGG": {
                    "status": "computed", "enrichment": 1.0,
                    "profile": [(10, 1.0), (-10, 2.0)],
                },
                "Z_TARGET": {
                    "status": "computed", "enrichment": 4.0,
                    "profile": [(10, 4.0), (-10, 2.0)],
                },
            },
            {
                "Z_TARGET": {
                    "status": "computed", "best_motif_id": "MA0139.1 <best>",
                    "best_adjusted_p_value": 0.001, "ame_status": "computed",
                },
            },
            {
                "A_IGG": [{"motif_id": "must-not-export"}],
                "Z_TARGET": [
                    {
                        "motif_id": f"M{rank:02d} & <motif>", "motif_alt_id": "CTCF",
                        "adjusted_p_value": rank / 1000, "p_value": rank / 100,
                        "effect": 2.0, "positive_sequences": 4, "rank": rank,
                    }
                    for rank in range(11, 0, -1)
                ],
            },
            annotation_status="computed_bed",
            insert_sizes={
                "A_IGG": [{"insert_size": 120, "pair_count": 3}],
                "Z_TARGET": [{"insert_size": 150, "pair_count": 5}],
            },
            peak_widths={
                "Z_TARGET": [{"width": 200, "peak_count": 2}],
            },
            fragments_per_peak={
                "Z_TARGET": [
                    {
                        "fragment_count": 4, "peak_count": 2,
                        "chrom": "chr1", "peak_name": "must-not-export",
                    },
                ],
            },
            motif_analysis_status="computed",
        )
        data["samples_by_id"]["Z_TARGET"]["warnings"].append(
            {"message": "Z_TARGET: NA warning & <visible>"}
        )
        data["warnings"] = [
            warning for sample in data["samples"] for warning in sample["warnings"]
        ]
        return data

    def test_representative_fixture_renders_twelve_samples_and_all_approved_panels(self):
        """The visual-QA fixture must exercise the complete offline dashboard."""
        workspace = self.make_workspace()
        output = workspace / "qc-dashboard-visual.html"
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tests" / "render_qc_dashboard_fixture.py"),
                str(output),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(output.is_file())
        html_document = output.read_text(encoding="utf-8")
        self.assertGreater(len(html_document), 10_000)
        self.assertTrue(html_document.startswith("<!doctype html>"))
        self.assertNotRegex(
            html_document,
            r'''(?:src|href)=["']https?://|url\(\s*["']?https?://''',
        )
        for sample_id in (
            "701_IgG", "702_IgG", "703_IgG",
            "701_CTCF", "702_CTCF", "703_CTCF",
            "704_GATA1", "705_GATA1", "706_GATA1",
            "704_RUNX1", "705_RUNX1", "706_RUNX1",
        ):
            with self.subTest(sample_id=sample_id):
                self.assertIn(f">{sample_id}<", html_document)
        for section_id in (
            "demultiplexing", "alignment", "insert-size-distribution",
            "peaks-frip", "peak-width-distribution", "tss-enrichment",
            "motif-enrichment",
        ):
            with self.subTest(section_id=section_id):
                self.assertIn(f'<section id="{section_id}">', html_document)
        for panel_title in (
            "Assigned read pairs", "Barcode balance within library",
            "Mapped reads", "Usable fragments after filtering",
            "PCR duplication", "End-to-end usable yield",
            "Peak count", "Fraction of reads in peaks",
            "Total bases covered by peaks", "Peak width median and range",
            "Peak count vs usable fragments", "Fragments per peak",
            "Insert-size distribution", "Peak-width distribution",
            "TSS enrichment score", "TSS profiles",
            "Motif enrichment heatmap",
        ):
            with self.subTest(panel_title=panel_title):
                self.assertIn(panel_title, html_document)

    def test_dashboard_derived_rows_use_approved_units_and_denominators(self):
        """A wrong denominator or silent zero would misstate library yield."""
        rows = qc.derive_dashboard_rows([
            {
                "sample_id": "S1",
                "sample_assigned_reads": 50,
                "sample_assignment_fraction": 0.625,
                "mapped_percent": 82.9,
                "mapq_filtered_fragments": 38,
                "duplicate_percent": 12,
                "peak_count": 4_321,
                "frip": 0.25,
                "total_covered_bases": 1_234_567,
            },
            {
                "sample_id": "ZERO",
                "sample_assigned_reads": 0,
                "mapq_filtered_fragments": 0,
            },
            {
                "sample_id": "MISSING",
                "sample_assigned_reads": None,
                "mapq_filtered_fragments": 10,
            },
        ])

        self.assertEqual(rows[0]["assigned_read_pairs_millions"], 0.00005)
        self.assertEqual(rows[0]["barcode_balance_percent"], 62.5)
        self.assertEqual(rows[0]["usable_fragments_millions"], 0.000038)
        self.assertEqual(rows[0]["end_to_end_yield_percent"], 76.0)
        self.assertEqual(rows[0]["peak_count_thousands"], 4.321)
        self.assertEqual(rows[0]["frip_percent"], 25.0)
        self.assertEqual(rows[0]["covered_megabases"], 1.234567)
        self.assertIsNone(rows[1]["end_to_end_yield_percent"])
        self.assertIsNone(rows[2]["end_to_end_yield_percent"])

    def test_dashboard_renders_all_approved_small_multiple_panels(self):
        """Dropping a panel would remove one of the approved cohort QC views."""
        data = self.report_data()
        target = data["samples_by_id"]["Z_TARGET"]
        control = data["samples_by_id"]["A_IGG"]
        target["library"].update({
            "mapped_percent": 82.9,
            "mapq_filtered_fragments": 38,
            "duplicate_percent": 12,
        })
        control["library"].update({
            "mapped_percent": 75,
            "mapq_filtered_fragments": 20,
            "duplicate_percent": 8,
        })
        target["peak"].update({
            "peak_count": 4_321,
            "total_covered_bases": 1_234_567,
            "width_min": 100,
            "width_q25": 125,
            "width_median": 170,
            "width_q75": 210,
            "width_max": 300,
            "frip": 0.25,
        })

        dashboard = qc.render_dashboard(data)

        for title in (
            "Assigned read pairs",
            "Barcode balance within library",
            "Mapped reads",
            "Usable fragments after filtering",
            "PCR duplication",
            "End-to-end usable yield",
            "Peak count",
            "Fraction of reads in peaks",
            "Total bases covered by peaks",
            "Peak width median and range",
            "Peak count vs usable fragments",
            "TSS enrichment score",
        ):
            self.assertIn(f"<h3>{title}</h3>", dashboard)
        self.assertIn('aria-label="Sequencing and alignment QC panels"', dashboard)
        self.assertIn('aria-label="Peak QC panels"', dashboard)
        self.assertIn('aria-label="TSS score panel"', dashboard)

    def test_technical_panels_include_controls_and_peak_panels_exclude_them(self):
        """IgG is technical QC data, but it has no biological peak call."""
        data = self.report_data()
        data["samples_by_id"]["A_IGG"]["library"].update({
            "mapped_percent": 75,
            "mapq_filtered_fragments": 20,
            "duplicate_percent": 8,
        })
        data["samples_by_id"]["Z_TARGET"]["peak"].update({
            "peak_count": 1,
            "total_covered_bases": 100,
            "width_min": 50,
            "width_q25": 60,
            "width_median": 70,
            "width_q75": 80,
            "width_max": 90,
            "frip": 0.2,
        })

        dashboard = qc.render_dashboard(data)
        sequencing = dashboard[
            dashboard.index('aria-label="Sequencing and alignment QC panels"'):
            dashboard.index('id="insert-size-distribution"')
        ]
        peaks = dashboard[
            dashboard.index('aria-label="Peak QC panels"'):
            dashboard.index('aria-label="Peaks and FRiP table"')
        ]

        self.assertIn('data-assay-target="IgG"', sequencing)
        self.assertNotIn('data-assay-target="IgG"', peaks)

    def test_serializers_use_stable_columns_and_json_null(self):
        """A rearranged or zero-filled summary would break downstream analysis."""
        workspace = self.make_workspace()
        data = self.report_data()

        qc.write_outputs(data, workspace)

        with (workspace / "qc_summary.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(list(rows[0]), qc.QC_SUMMARY_COLUMNS)
        self.assertTrue({
            "peak_width_min", "peak_width_q25", "peak_width_mean",
            "peak_width_median", "peak_width_q75", "peak_width_max",
        }.issubset(rows[0]))
        self.assertEqual(rows[0]["frip"], "")
        self.assertEqual([row["sample_id"] for row in rows], ["A_IGG", "Z_TARGET"])
        self.assertEqual(rows[0]["is_control"], "true")
        self.assertEqual(rows[1]["is_control"], "false")
        payload = json.loads((workspace / "qc_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(
            set(payload),
            {
                "schema_version", "generator_version", "annotation_status",
                "counts", "availability", "metric_definitions", "samples",
                "warnings",
            },
        )
        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(payload["generator_version"], "1.1.0")
        self.assertEqual(
            payload["counts"],
            {"samples": 2, "targets": 1, "controls": 1, "warnings": 3},
        )
        self.assertEqual(
            set(payload["availability"]),
            {
                "demultiplex", "library", "insert_size", "peak",
                "peak_width", "tss", "motif", "ame",
            },
        )
        self.assertEqual(
            set(payload["availability"]["tss"]),
            {
                "status", "computed", "skipped", "empty", "missing",
                "failed", "not_applicable",
            },
        )
        self.assertEqual(
            set(payload["samples"][0]),
            {
                "sample_id", "library_id", "input_group", "assay_target",
                "is_control", "control_id", "expected_motif", "sample_kind",
                "demultiplex", "library", "peak", "tss", "motif",
                "availability", "top_motifs", "warnings",
            },
        )
        self.assertIsNone(payload["samples"][0]["peak"]["frip"])
        self.assertIsNone(
            payload["samples"][0]["peak"]["fragments_per_peak_distribution"]
        )
        self.assertEqual(
            payload["metric_definitions"]["tss_enrichment"]["formula"],
            "center_bin_signal / mean(terminal_100bp_flanks)",
        )

    def test_json_whitelist_blocks_arbitrary_producer_columns(self):
        """Adding an internal producer column must not silently expand schema v1."""
        workspace = self.make_workspace()
        data = self.report_data()
        data["samples_by_id"]["Z_TARGET"]["library"]["future_private_metric"] = 9

        qc.write_qc_summary_json(data, workspace / "qc_summary.json")

        payload = json.loads((workspace / "qc_summary.json").read_text(encoding="utf-8"))
        target = next(row for row in payload["samples"] if row["sample_id"] == "Z_TARGET")
        self.assertNotIn("future_private_metric", target["library"])

    def test_distribution_inputs_change_public_json_and_dashboard(self):
        """Both staged histogram families must affect consumer-visible outputs."""
        workspace = self.make_workspace()
        data = self.report_data()
        data["samples_by_id"]["Z_TARGET"]["library"]["insert_size_distribution"] = [
            {"insert_size": 147, "pair_count": 23},
        ]
        data["samples_by_id"]["Z_TARGET"]["peak"]["width_distribution"] = [
            {"width": 321, "peak_count": 7},
        ]
        data["samples_by_id"]["Z_TARGET"]["availability"]["insert_size"] = {
            "status": "computed", "reason": None,
        }
        data["samples_by_id"]["Z_TARGET"]["availability"]["peak_width"] = {
            "status": "computed", "reason": None,
        }

        qc.write_outputs(data, workspace)

        payload = json.loads((workspace / "qc_summary.json").read_text(encoding="utf-8"))
        target = next(row for row in payload["samples"] if row["sample_id"] == "Z_TARGET")
        self.assertEqual(
            target["library"]["insert_size_distribution"],
            [{"insert_size": 147, "pair_count": 23}],
        )
        self.assertEqual(
            target["peak"]["width_distribution"],
            [{"width": 321, "peak_count": 7}],
        )
        self.assertEqual(
            target["peak"]["fragments_per_peak_distribution"],
            [{"fragment_count": 4, "peak_count": 2}],
        )
        self.assertNotIn(
            "chrom",
            target["peak"]["fragments_per_peak_distribution"][0],
        )
        self.assertNotIn(
            "peak_name",
            target["peak"]["fragments_per_peak_distribution"][0],
        )
        dashboard = (workspace / "qc_dashboard.html").read_text(encoding="utf-8")
        self.assertIn("Insert-size distribution", dashboard)
        self.assertIn("250 bp bins", dashboard)
        self.assertIn(">[0, 250)<", dashboard)
        self.assertIn(">23.0<", dashboard)
        self.assertIn("Peak-width distribution", dashboard)
        self.assertIn(">[250, 500)<", dashboard)
        self.assertIn(">7.00<", dashboard)

    def test_dashboard_bins_insert_and_peak_width_distributions_across_samples(self):
        data = self.report_data()
        data["samples_by_id"]["A_IGG"]["library"]["insert_size_distribution"] = [
            {"insert_size": 10, "pair_count": 1},
        ]
        data["samples_by_id"]["Z_TARGET"]["library"]["insert_size_distribution"] = [
            {"insert_size": 260, "pair_count": 1},
            {"insert_size": 510, "pair_count": 2},
        ]
        data["samples_by_id"]["Z_TARGET"]["peak"]["width_distribution"] = [
            {"width": 20, "peak_count": 1},
            {"width": 260, "peak_count": 3},
        ]

        dashboard = qc.render_dashboard(data)

        insert_section = dashboard[
            dashboard.index('id="insert-size-distribution"'):
            dashboard.index('id="peaks-frip"')
        ]
        self.assertIn("Percent of read pairs", insert_section)
        self.assertIn(">[0, 250)<", insert_section)
        self.assertIn(">[250, 500)<", insert_section)
        self.assertIn(">[500, 750)<", insert_section)
        self.assertIn(">100<", insert_section)
        self.assertAlmostEqual(
            sum(
                float(value)
                for value in re.findall(
                    r'data-sample-id="Z_TARGET"[^>]*data-percent="([0-9.]+)"',
                    insert_section,
                )
            ),
            100.0,
        )

        width_section = dashboard[
            dashboard.index('id="peak-width-distribution"'):
            dashboard.index('id="tss-enrichment"')
        ]
        self.assertIn("Percent of peaks", width_section)
        self.assertIn(">[0, 250)<", width_section)
        self.assertIn(">[250, 500)<", width_section)

    def test_dashboard_renders_fragments_per_peak_ecdf_and_target_colored_tss_profiles(self):
        data = self.report_data()
        data["samples_by_id"]["Z_TARGET"]["peak"][
            "fragments_per_peak_distribution"
        ] = [
            {"fragment_count": 0, "peak_count": 2},
            {"fragment_count": 1, "peak_count": 1},
            {"fragment_count": 10, "peak_count": 1},
        ]

        dashboard = qc.render_dashboard(data)

        peak_section = dashboard[
            dashboard.index('id="peaks-frip"'):
            dashboard.index('id="peak-width-distribution"')
        ]
        self.assertIn("<h3>Fragments per peak</h3>", peak_section)
        self.assertIn("Cumulative percent of peaks", peak_section)
        self.assertIn("50.0% of peaks have zero fragments", peak_section)
        self.assertIn('data-zero-origin="true"', peak_section)

        tss_section = dashboard[
            dashboard.index('id="tss-enrichment"'):
            dashboard.index('id="motif-enrichment"')
        ]
        self.assertIn("Mean coverage (RPKM)", tss_section)
        self.assertIn('class="zero-reference"', tss_section)
        self.assertIn('data-assay-target="CTCF"', tss_section)
        self.assertIn('stroke="#2F78D1"', tss_section)
        self.assertIn('class="endpoint-label"', tss_section)

    def test_dashboard_renders_target_only_heatmap_and_displayed_cell_table(self):
        data = self.report_data()
        target = data["samples_by_id"]["Z_TARGET"]
        records = [
            {
                "motif_id": f"MA_OTHER_{rank:02d}",
                "motif_alt_id": f"OTHER{rank:02d}",
                "rank": rank,
                "adjusted_p_value": rank / 1000,
                "p_value": rank / 100,
                "effect": 2.0,
                "positive_sequences": 4,
            }
            for rank in range(1, 12)
        ]
        target["ame_motifs"] = records
        target["top_motifs"] = records[:10]
        data["samples_by_id"]["A_IGG"]["ame_motifs"] = [{
            "motif_id": "CONTROL_ONLY",
            "motif_alt_id": "CONTROL",
            "rank": 1,
            "adjusted_p_value": 0,
        }]

        dashboard = qc.render_dashboard(data)
        motif_section = dashboard[
            dashboard.index('id="motif-enrichment"'):
            dashboard.index('id="warnings"')
        ]

        self.assertIn('aria-label="Motif enrichment heatmap"', motif_section)
        self.assertIn("MA_OTHER_11", motif_section)
        self.assertNotIn("CONTROL_ONLY", motif_section)
        self.assertNotIn(">A_IGG<", motif_section)
        compact_table = re.search(
            r'aria-label="Displayed motif cells table".*?</table>',
            motif_section,
            re.DOTALL,
        )
        self.assertIsNotNone(compact_table)
        self.assertEqual(compact_table.group().count("<tr>"), 12)
        for heading in (
            "Sample", "Motif", "Alternate ID", "Rank",
            "Adjusted p-value", "Transformed significance", "Cognate",
        ):
            self.assertIn(heading, compact_table.group())

    def test_eleventh_heatmap_motif_does_not_expand_top_motifs_export(self):
        workspace = self.make_workspace()
        data = self.report_data()
        records = [
            {
                "motif_id": f"MA_OTHER_{rank:02d}",
                "motif_alt_id": f"OTHER{rank:02d}",
                "rank": rank,
                "adjusted_p_value": rank / 1000,
                "p_value": rank / 100,
                "effect": 2.0,
                "positive_sequences": 4,
            }
            for rank in range(1, 12)
        ]
        data["samples_by_id"]["Z_TARGET"]["ame_motifs"] = records
        data["samples_by_id"]["Z_TARGET"]["top_motifs"] = records

        dashboard = qc.render_dashboard(data)
        qc.write_top_motifs_tsv(data, workspace / "top_motifs.tsv")

        motif_section = dashboard[
            dashboard.index('id="motif-enrichment"'):
            dashboard.index('id="warnings"')
        ]
        self.assertIn("MA_OTHER_11", motif_section)
        with (workspace / "top_motifs.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            exported = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(exported), 10)
        self.assertNotIn("MA_OTHER_11", {row["motif_id"] for row in exported})

    def test_tidy_exports_exclude_controls_limit_motifs_and_order_rows(self):
        """Controls must not leak into motif biology and sortable exports stay stable."""
        workspace = self.make_workspace()
        data = self.report_data()
        data["samples_by_id"]["Z_TARGET"]["top_motifs"] = [
            {"motif_id": "RANK_2_LOW_P", "adjusted_p_value": 0.000001, "rank": 2},
            {"motif_id": "B_RANK_1", "adjusted_p_value": 0.9, "rank": 1},
            {"motif_id": "A_RANK_1", "adjusted_p_value": 0.8, "rank": 1},
        ]

        qc.write_outputs(data, workspace)

        with (workspace / "top_motifs.tsv").open(encoding="utf-8", newline="") as handle:
            motif_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(motif_rows), 3)
        self.assertEqual({row["sample_id"] for row in motif_rows}, {"Z_TARGET"})
        self.assertEqual(
            [row["motif_id"] for row in motif_rows],
            ["A_RANK_1", "B_RANK_1", "RANK_2_LOW_P"],
        )
        with (workspace / "tss_profiles.tsv").open(encoding="utf-8", newline="") as handle:
            profile_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(list(profile_rows[0]), ["sample_id", "position_bp", "signal"])
        self.assertEqual(
            [(row["sample_id"], int(row["position_bp"])) for row in profile_rows],
            [("A_IGG", -10), ("A_IGG", 10), ("Z_TARGET", -10), ("Z_TARGET", 10)],
        )

    def test_dashboard_is_offline_semantic_descriptive_and_escaped(self):
        """Unescaped input or a hidden classification changes the report contract."""
        html = qc.render_dashboard(self.report_data())

        for section_id in (
            "run-overview", "demultiplexing", "alignment", "peaks-frip",
            "tss-enrichment", "motif-enrichment",
        ):
            self.assertIn(f'id="{section_id}"', html)
        self.assertIn("<svg", html)
        self.assertNotIn("https://", html)
        self.assertNotIn("http://", html)
        self.assertNotIn("overall pass", html.lower())
        self.assertNotIn("overall fail", html.lower())
        self.assertNotIn('class="qc-pass"', html)
        self.assertNotIn('class="qc-fail"', html)
        for visible_text in (
            "NA", "Targets (solid)", "IgG controls (outlined)",
            "Bar charts — Targets (solid)",
            "Line charts use the per-series swatches",
            "CTCF &amp; &lt;target&gt;", "M01 &amp; &lt;motif&gt;",
            "NA warning &amp; &lt;visible&gt;",
        ):
            self.assertIn(visible_text, html)
        self.assertNotIn("CTCF & <target>", html)
        self.assertNotIn("M01 & <motif>", html)
        control_trace = re.search(r'<polyline[^>]*data-sample-kind="control"[^>]*>', html)
        target_trace = re.search(r'<polyline[^>]*data-sample-kind="target"[^>]*>', html)
        self.assertIsNotNone(control_trace)
        self.assertIsNotNone(target_trace)
        self.assertIn('data-marker="square"', control_trace.group())
        self.assertIn('data-marker="circle"', target_trace.group())

    def test_dashboard_renders_approved_counts_scalars_control_peak_na_and_motif_status(self):
        """Scalar fields already joined into the model must not disappear in HTML."""
        data = self.report_data()
        target = data["samples_by_id"]["Z_TARGET"]
        target["demultiplex"].update({
            "total_read_pairs": 101,
            "assigned_read_pairs": 80,
            "unassigned_read_pairs": 16,
            "sample_assigned_reads": 50,
        })
        target["library"].update({
            "raw_total_reads": 111,
            "mapq_filtered_reads": 77,
            "mapq_filtered_fragments": 38,
            "mapq_filtered_fraction": 0.693,
            "markdup_examined_reads": 100,
            "duplicate_total": 12,
            "duplicate_percent": 12,
            "mitochondrial_percent": 1.5,
            "estimated_library_size": 900,
            "insert_size_total_pairs": 42,
            "insert_size_min": 90,
            "insert_size_q25": 110,
            "insert_size_mean": 140,
            "insert_size_median": 135,
            "insert_size_q75": 160,
            "insert_size_max": 220,
        })
        target["peak"].update({
            "peak_count": 4,
            "total_covered_bases": 1234,
            "width_min": 100,
            "width_q25": 125,
            "width_mean": 175,
            "width_median": 170,
            "width_q75": 210,
            "width_max": 300,
            "total_fragments": 40,
            "fragments_in_peaks": 10,
            "frip": 0.25,
        })

        dashboard = qc.render_dashboard(data)

        for text in (
            "2 samples", "1 target", "1 control", "Input-family availability",
            "Total read pairs", "Assigned read pairs", "Sample assigned reads",
            "Raw reads", "MAPQ reads", "MAPQ fragments", "MAPQ fraction",
            "Examined reads", "Duplicate reads", "Duplicate (%)",
            "Mitochondrial (%)", "Estimated library size", "Insert pairs",
            "Covered bases", "Width min", "Width Q25", "Width mean",
            "Width median", "Width Q75", "Width max", "Usable fragments",
            "Fragments in peaks", "Best adjusted significance", "AME status",
        ):
            self.assertIn(text, dashboard)
        peak_table = re.search(
            r'aria-label="Peaks and FRiP table".*?</table>',
            dashboard,
            re.DOTALL,
        )
        self.assertIsNotNone(peak_table)
        self.assertRegex(peak_table.group(), r"A_IGG.*?NA.*?NA")
        self.assertIn("0.001", dashboard)
        self.assertIn("computed", dashboard)

    def test_run_overview_retains_library_group_control_and_motif_metadata(self):
        """The overview must show the matched-control context needed to compare samples."""
        dashboard = qc.render_dashboard(self.report_data())
        overview = re.search(
            r'aria-label="Run overview table".*?</table>',
            dashboard,
            re.DOTALL,
        )

        self.assertIsNotNone(overview)
        for text in (
            "Physical library", "Input group", "Matched control",
            "Expected motif", "A_IGG", "CTCF &amp; &lt;target&gt;",
        ):
            self.assertIn(text, overview.group())

    def test_dashboard_tables_are_scrollable_without_page_width_overflow(self):
        """Narrow viewports need per-table scrolling instead of a clipped page."""
        html = qc.render_dashboard(self.report_data())

        self.assertGreater(html.count("<table>"), 0)
        self.assertEqual(
            html.count("<table>"),
            html.count(
                '<div class="table-scroll" role="region" '
                'aria-label="'
            ),
        )
        self.assertIn(
            ".table-scroll{max-width:100%;overflow-x:auto}",
            html,
        )

    def test_dashboard_rounds_only_html_and_collapses_detailed_tables(self):
        """Presentation compaction must not alter exact identifiers or raw exports."""
        workspace = self.make_workspace()
        data = self.report_data()
        target = data["samples_by_id"]["Z_TARGET"]
        target["demultiplex"].update({
            "assigned_read_pairs": 17_432_100,
            "assigned_fraction": 0.012345,
        })
        target["library"]["insert_size_distribution"] = [
            {"insert_size": 260, "pair_count": 17_432_100},
        ]
        target["availability"]["insert_size"] = {
            "status": "computed", "reason": None,
        }

        qc.write_outputs(data, workspace)

        dashboard = (workspace / "qc_dashboard.html").read_text(encoding="utf-8")
        self.assertIn(">17.4M<", dashboard)
        self.assertIn(">0.0123<", dashboard)
        self.assertIn("<details", dashboard)
        self.assertIn(
            "<summary>View binned insert-size data</summary>",
            dashboard,
        )
        self.assertNotIn("<details open", dashboard)
        self.assertIn(">Z_TARGET<", dashboard)
        self.assertIn(">1<", dashboard)
        self.assertIn(">[250, 500)<", dashboard)

        with (workspace / "qc_summary.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        raw_target = next(row for row in rows if row["sample_id"] == "Z_TARGET")
        self.assertEqual(raw_target["assigned_read_pairs"], "17432100")
        self.assertEqual(raw_target["assigned_fraction"], "0.012345")

        payload = json.loads(
            (workspace / "qc_summary.json").read_text(encoding="utf-8")
        )
        json_target = next(
            row for row in payload["samples"]
            if row["sample_id"] == "Z_TARGET"
        )
        self.assertEqual(
            json_target["demultiplex"]["assigned_read_pairs"],
            17_432_100,
        )
        self.assertEqual(
            json_target["demultiplex"]["assigned_fraction"],
            0.012345,
        )

    def test_dashboard_version_and_responsive_layout_metadata(self):
        """The redesigned HTML and producer record must identify version 1.1.0."""
        dashboard_module = (
            Path(__file__).parents[2] / "modules/local/qc_dashboard.nf"
        ).read_text(encoding="utf-8")
        rendered = qc.render_dashboard(self.report_data())

        self.assertEqual(qc.SCHEMA_VERSION, 1)
        self.assertEqual(qc.GENERATOR_VERSION, "1.1.0")
        self.assertIn("qc_dashboard.py: 1.1.0", dashboard_module)
        self.assertIn(".panel-grid{", rendered)
        self.assertIn("@media (max-width:", rendered)
        self.assertIn("@media print{", rendered)

    def test_dashboard_tables_have_unique_context_specific_accessible_names(self):
        """Assistive technology must distinguish every scrollable data region."""
        html = qc.render_dashboard(self.report_data())

        names = re.findall(
            r'<div class="table-scroll" role="region" '
            r'aria-label="([^"]+)" tabindex="0">',
            html,
        )
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual(
            set(names),
            {
                "Run overview table",
                "Input-family availability table",
                "Demultiplexing metrics table",
                "Alignment and library QC table",
                "Insert-size distribution table",
                "Peaks and FRiP table",
                "Peak-width distribution table",
                "TSS enrichment table",
                "Expected motif table",
                "Top AME motifs table",
                "Warnings table",
            },
        )

    def test_bar_chart_scales_to_one_hundred_long_labels_without_overlap(self):
        """A fixed plot width must not collapse realistic cohorts into overlapping bars."""
        chart = qc.render_bar_chart(
            "Mapped reads",
            [
                {
                    "sample_id": f"long_sample_identifier_{index:03d}_replicate_alpha",
                    "mapped_percent": index + 1,
                    "is_control": index % 10 == 0,
                }
                for index in range(100)
            ],
            value_key="mapped_percent",
        )

        viewbox = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', chart)
        self.assertIsNotNone(viewbox)
        self.assertGreater(float(viewbox.group(1)), 640)
        self.assertGreater(float(viewbox.group(2)), 400)
        bars = [
            (float(x), float(width))
            for x, width in re.findall(
                r'<rect class="bar" x="([0-9.]+)"[^>]*width="([0-9.]+)"',
                chart,
            )
        ]
        self.assertEqual(len(bars), 100)
        self.assertTrue(all(
            left_x + left_width <= right_x
            for (left_x, left_width), (right_x, _) in zip(bars, bars[1:])
        ))
        self.assertIn(
            '<div class="chart-scroll" role="region" '
            'aria-label="Mapped reads chart" tabindex="0">',
            chart,
        )
        self.assertIn("long_sample_identifier_099_replicate_alpha", chart)

    def test_line_chart_uses_fixed_plot_and_unique_non_color_mapping_for_many_targets(self):
        """More than five TSS profiles need a textual/style key independent of color."""
        chart = qc.render_line_chart(
            "TSS profiles",
            {
                f"long_target_sample_{index:02d}": {
                    "profile": [(-10, 1), (10, index + 1)],
                    "is_control": index == 0,
                }
                for index in range(12)
            },
        )

        viewbox = re.search(r'viewBox="0 0 640 ([0-9.]+)"', chart)
        self.assertIsNotNone(viewbox)
        self.assertLessEqual(float(viewbox.group(1)), 320)
        traces = re.findall(
            r'<polyline[^>]*data-series-key="([^"]+)"[^>]*'
            r'data-sample-kind="([^"]+)"[^>]*'
            r'data-marker="([^"]+)"[^>]*'
            r'stroke-dasharray="([^"]+)"',
            chart,
        )
        self.assertEqual(len(traces), 12)
        self.assertEqual(len({key for key, _, _, _ in traces}), 12)
        self.assertEqual(len({dash for _, _, _, dash in traces}), 12)
        self.assertEqual(traces[0][1:3], ("control", "square"))
        self.assertTrue(
            all(marker == "circle" for _, kind, marker, _ in traces if kind == "target")
        )
        legend_pairs = re.findall(
            r'<span class="series-key">([^<]+)</span>'
            r'<span class="series-label">([^<]+)</span>',
            chart,
        )
        self.assertEqual(len(legend_pairs), 12)
        self.assertEqual(
            {label for _, label in legend_pairs},
            {f"long_target_sample_{index:02d}" for index in range(12)},
        )
        self.assertEqual(chart.count('class="series-swatch"'), 12)

    def test_distribution_chart_uses_a_fixed_numeric_axis_for_large_histograms(self):
        """Distribution geometry must scale by series, not by every sample-bin pair."""
        rows = [
            {
                "sample_id": f"sample_{sample_index:03d}",
                "is_control": sample_index % 10 == 0,
                "insert_size": bin_index,
                "pair_count": (sample_index + bin_index) % 101,
            }
            for sample_index in range(100)
            for bin_index in range(1, 501)
        ]

        chart = qc.render_distribution_chart(
            "Insert-size distribution",
            rows,
            x_key="insert_size",
            value_key="pair_count",
            x_axis_label="Insert size (bp)",
            y_axis_label="Read pairs",
        )

        self.assertIn('viewBox="0 0 640 280"', chart)
        traces = re.findall(
            r'<polyline[^>]*data-series-key="([^"]+)"',
            chart,
        )
        self.assertEqual(len(traces), 100)
        self.assertEqual(len(set(traces)), 100)
        self.assertEqual(chart.count('class="series-swatch"'), 100)
        self.assertEqual(chart.count('class="axis-tick-label"'), 6)
        self.assertLess(len(chart), 2_000_000)

    def test_distribution_chart_centers_a_single_bin_with_one_x_tick(self):
        """A degenerate numeric range must not stack labels at the left edge."""
        chart = qc.render_distribution_chart(
            "Single-bin distribution",
            [{
                "sample_id": "sample_001",
                "is_control": False,
                "insert_size": 150,
                "pair_count": 12,
            }],
            x_key="insert_size",
            value_key="pair_count",
            x_axis_label="Insert size (bp)",
            y_axis_label="Read pairs",
        )

        self.assertEqual(chart.count('class="axis-tick-label"'), 4)
        self.assertRegex(chart, r'<polyline[^>]*points="338\.0,[0-9.]+?"')

    def test_cli_publishes_a_complete_output_set_and_reports_input_errors(self):
        """A partial report must never be published after a bad command invocation."""
        workspace = self.make_workspace()
        metadata = workspace / "sample_metadata.json"
        metadata.write_text(json.dumps([
            {
                "sample_id": "S1", "library_id": "L1", "input_group": "25K",
                "assay_target": "CTCF", "is_control": False, "control_id": "I1",
                "expected_motif": "CTCF",
            },
            {
                "sample_id": "I1", "library_id": "L1", "input_group": "25K",
                "assay_target": "IgG", "is_control": True, "control_id": None,
                "expected_motif": None,
            },
        ]), encoding="utf-8")
        directories = {}
        for name in ("demux", "library", "insert", "peak", "tss", "motif", "ame"):
            directories[name] = workspace / name
            directories[name].mkdir()
        (directories["demux"] / "L1.metrics.json").write_text(json.dumps({
            "total_reads": 4, "assigned_reads": 4, "ambiguous_reads": 0,
            "unassigned_reads": 0, "assignment_counts": {"I1": 1, "S1": 3},
        }), encoding="utf-8")
        (directories["library"] / "S1.library_qc.tsv").write_text(
            library_metrics_tsv("S1", mapped_percent=95), encoding="utf-8"
        )
        (directories["library"] / "I1.library_qc.tsv").write_text(
            library_metrics_tsv("I1", mapped_percent=90), encoding="utf-8"
        )
        for sample_id, insert_size, pair_count in (("S1", 147, 3), ("I1", 121, 1)):
            (directories["insert"] / f"{sample_id}.insert_size_distribution.tsv").write_text(
                "sample_id\tinsert_size\tpair_count\n"
                f"{sample_id}\t{insert_size}\t{pair_count}\n",
                encoding="utf-8",
            )
            shutil.copyfile(
                GOLDEN_DEEPTOOLS_PROFILE,
                directories["tss"] / f"{sample_id}.tss_profile.tsv",
            )
            (directories["tss"] / f"{sample_id}.tss_status.tsv").write_text(
                "sample_id\tannotation_mode\tstatus\n"
                f"{sample_id}\tbed\tcomputed\n",
                encoding="utf-8",
            )
        (directories["peak"] / "S1.peak_qc.tsv").write_text(
            peak_metrics_tsv("S1"),
            encoding="utf-8",
        )
        (directories["peak"] / "S1.peak_qc.width_histogram.tsv").write_text(
            "width\tpeak_count\n321\t1\n",
            encoding="utf-8",
        )
        outdir = workspace / "dashboard"
        arguments = [
            "--metadata", str(metadata), "--demux-dir", str(directories["demux"]),
            "--library-dir", str(directories["library"]), "--insert-dir", str(directories["insert"]),
            "--peak-dir", str(directories["peak"]), "--tss-dir", str(directories["tss"]),
            "--motif-dir", str(directories["motif"]), "--ame-dir", str(directories["ame"]),
            "--annotation-status", "computed_bed",
            "--motif-analysis-status", "skipped_no_database",
            "--outdir", str(outdir),
        ]

        self.assertEqual(qc.main(arguments), 0)
        self.assertEqual(
            {path.name for path in outdir.iterdir()},
            {"qc_dashboard.html", "qc_summary.tsv", "qc_summary.json", "top_motifs.tsv", "tss_profiles.tsv"},
        )
        with (outdir / "qc_summary.tsv").open(encoding="utf-8", newline="") as handle:
            rows = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(rows["S1"]["tss_enrichment"], "6")
        payload = json.loads((outdir / "qc_summary.json").read_text(encoding="utf-8"))
        target = next(row for row in payload["samples"] if row["sample_id"] == "S1")
        self.assertEqual(
            target["library"]["insert_size_distribution"],
            [{"insert_size": 147, "pair_count": 3}],
        )
        self.assertEqual(
            target["peak"]["width_distribution"],
            [{"width": 321, "peak_count": 1}],
        )
        error = io.StringIO()
        with contextlib.redirect_stderr(error):
            result = qc.main(["--metadata", str(workspace / "missing.json"), *arguments[2:]])
        self.assertEqual(result, 2)
        self.assertIn("qc_dashboard.py: error:", error.getvalue())
        self.assertTrue((outdir / "qc_dashboard.html").is_file())

    def test_atomic_bundle_publication_rolls_back_after_each_directory_replace(self):
        """Injected publication faults must leave the complete prior five-file set."""
        workspace = self.make_workspace()
        names = (
            "qc_dashboard.html", "qc_summary.tsv", "qc_summary.json",
            "top_motifs.tsv", "tss_profiles.tsv",
        )
        data = {
            "schema_version": 1,
            "annotation_status": "skipped_no_annotation",
            "samples": [],
            "warnings": [],
        }
        for fail_after in (1, 2):
            with self.subTest(fail_after=fail_after):
                outdir = workspace / f"dashboard-{fail_after}"
                outdir.mkdir()
                for name in names:
                    (outdir / name).write_bytes(f"OLD::{name}\n".encode())
                calls = 0

                def replace_then_fail(source, destination):
                    nonlocal calls
                    os.replace(source, destination)
                    calls += 1
                    if calls == fail_after:
                        raise OSError(f"injected fault after replace {fail_after}")

                with mock.patch.object(
                    qc, "_replace_path", side_effect=replace_then_fail, create=True
                ):
                    with self.assertRaisesRegex(
                        OSError, f"injected fault after replace {fail_after}"
                    ):
                        qc.write_outputs_atomically(data, outdir)

                self.assertEqual(
                    {path.name: path.read_bytes() for path in outdir.iterdir()},
                    {name: f"OLD::{name}\n".encode() for name in names},
                )
                self.assertEqual(
                    list(workspace.glob(f".{outdir.name}.*")),
                    [],
                )

    def test_atomic_bundle_refuses_to_replace_a_nondedicated_directory(self):
        """Directory publication must never delete an unrelated caller-owned file."""
        workspace = self.make_workspace()
        outdir = workspace / "shared"
        outdir.mkdir()
        sentinel = outdir / "keep-me.txt"
        sentinel.write_bytes(b"caller-owned\n")
        data = {
            "schema_version": 1,
            "annotation_status": "skipped_no_annotation",
            "samples": [],
            "warnings": [],
        }

        with self.assertRaisesRegex(
            qc.DashboardInputError,
            "dedicated dashboard output directory",
        ):
            qc.write_outputs_atomically(data, outdir)
        with self.assertRaisesRegex(
            qc.DashboardInputError,
            "dedicated dashboard output directory",
        ):
            qc.write_outputs_atomically(data, Path.cwd())

        self.assertEqual(
            {path.name: path.read_bytes() for path in outdir.iterdir()},
            {"keep-me.txt": b"caller-owned\n"},
        )
        self.assertEqual(list(workspace.glob(".shared.*")), [])


if __name__ == "__main__":
    unittest.main()
