import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

import qc_dashboard as qc


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
            "assignment_counts": {"S1": 30, "S2": 50},
        }), encoding="utf-8")

        rows = qc.read_demultiplex_metrics([path])

        self.assertEqual(rows["L1"]["assigned_fraction"], 0.8)
        self.assertEqual(rows["L1"]["assignment_counts"]["S1"], 30)

    def test_read_demultiplex_metrics_rejects_duplicate_library_ids(self):
        """A library must contribute exactly one unambiguous count record."""
        workspace = self.make_workspace()
        paths = []
        for index in range(2):
            path = workspace / f"demultiplex-{index}.json"
            path.write_text(json.dumps({
                "library_id": "L1", "total_reads": 1, "assigned_reads": 1,
                "ambiguous_reads": 0, "unassigned_reads": 0, "assignment_counts": {"S1": 1},
            }), encoding="utf-8")
            paths.append(path)

        with self.assertRaisesRegex(qc.DashboardInputError, "duplicate library_id L1"):
            qc.read_demultiplex_metrics(paths)

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
                    qc.read_demultiplex_metrics([path])

    def test_read_library_metrics_rejects_duplicate_sample_ids(self):
        """A sample may have only one library-QC summary row."""
        workspace = self.make_workspace()
        paths = []
        for index in range(2):
            path = workspace / f"library-{index}.tsv"
            path.write_text(
                "sample_id\traw_total_reads\tmapped_percent\nS1\t100\t95\n",
                encoding="utf-8",
            )
            paths.append(path)

        with self.assertRaisesRegex(qc.DashboardInputError, "duplicate sample_id S1"):
            qc.read_library_metrics(paths)

    def test_read_peak_metrics_parses_frip_and_rejects_duplicate_metric_names(self):
        """Metric/value peak tables require one value per metric name."""
        workspace = self.make_workspace()
        valid = workspace / "S1.peak_qc.tsv"
        valid.write_text(
            "metric\tvalue\nsample_id\tS1\npeak_count\t12\nfrip\t0.42\n",
            encoding="utf-8",
        )
        self.assertEqual(qc.read_peak_metrics([valid])["S1"]["frip"], 0.42)

        duplicate = workspace / "duplicate.peak_qc.tsv"
        duplicate.write_text(
            "metric\tvalue\nsample_id\tS2\nfrip\t0.42\nfrip\t0.43\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(qc.DashboardInputError, "duplicate peak metric frip"):
            qc.read_peak_metrics([duplicate])


class TssAndAmeTests(unittest.TestCase):
    def make_workspace(self):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        return Path(temporary_directory.name)

    def test_tss_enrichment_uses_center_and_terminal_100bp(self):
        """The scalar score must use the center and both terminal 100-bp flanks."""
        workspace = self.make_workspace()
        values = [2.0] * 600
        values[300] = 12.0
        path = workspace / "S1.tss_profile.tsv"
        path.write_text("sample\tgroup\t" + "\t".join(map(str, values)) + "\n", encoding="utf-8")

        profile = qc.read_tss_profile(path)

        self.assertEqual(profile[0], (-3000, 2.0))
        self.assertEqual(profile[300], (0, 12.0))
        self.assertEqual(qc.calculate_tss_enrichment(profile), 6.0)

    def test_tss_zero_flank_returns_missing(self):
        """A zero baseline has no defined center-to-flank enrichment ratio."""
        profile = [(index * 10 - 3000, 0.0) for index in range(600)]

        self.assertIsNone(qc.calculate_tss_enrichment(profile))

    def test_read_tss_profile_rejects_wrong_bin_count_and_nonfinite_values(self):
        """A profile with an incompatible schema cannot be joined to the report."""
        workspace = self.make_workspace()
        wrong_count = workspace / "wrong-count.tsv"
        wrong_count.write_text("sample\t" + "\t".join(["1"] * 599) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(qc.DashboardInputError, "600 numeric bins"):
            qc.read_tss_profile(wrong_count)

        nonfinite = workspace / "nonfinite.tsv"
        nonfinite.write_text("sample\t" + "\t".join(["1"] * 599 + ["nan"]) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(qc.DashboardInputError, "must be finite"):
            qc.read_tss_profile(nonfinite)

    def test_read_tss_profile_rejects_missing_or_multiple_profiles(self):
        """The dashboard needs exactly one aggregate profile row per input file."""
        workspace = self.make_workspace()
        missing = workspace / "missing.tsv"
        missing.write_text("# no aggregate profile\n\n", encoding="utf-8")
        with self.assertRaisesRegex(qc.DashboardInputError, "no TSS profile row"):
            qc.read_tss_profile(missing)

        multiple = workspace / "multiple.tsv"
        row = "sample\t" + "\t".join(["1"] * 600) + "\n"
        multiple.write_text(row + row, encoding="utf-8")
        with self.assertRaisesRegex(qc.DashboardInputError, "multiple compatible"):
            qc.read_tss_profile(multiple)

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


if __name__ == "__main__":
    unittest.main()
