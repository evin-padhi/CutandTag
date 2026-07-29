import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

import peak_qc
from peak_qc import (
    Fragment,
    Peak,
    PeakQcError,
    calculate_frip,
    read_peaks,
    summarize_peaks,
    write_qc_outputs,
)


class PeakQcUnitTests(unittest.TestCase):
    def setUp(self):
        self.peaks = [
            ("chr1", 10, 20, "peak_a"),
            ("chr1", 15, 30, "peak_b"),
            ("chr2", 0, 5, "peak_c"),
        ]
        self.fragments = [
            ("chr1", 5, 10, "chr1", 20, 25, "fragment_one"),
            ("chr1", 19, 21, "chr1", 21, 22, "fragment_two"),
            ("chr1", 30, 31, "chr1", 31, 32, "fragment_three"),
            ("chr2", 4, 6, "chr2", 6, 7, "fragment_four"),
        ]

    def test_peak_summary_reports_union_coverage_and_width_quantiles(self):
        summary = summarize_peaks(self.peaks)

        self.assertEqual(summary["peak_count"], 3)
        self.assertEqual(summary["total_covered_bases"], 25)
        self.assertEqual(summary["peak_width_min"], 5)
        self.assertEqual(summary["peak_width_mean"], 10)
        self.assertEqual(summary["peak_width_median"], 10)
        self.assertEqual(summary["peak_width_max"], 15)
        self.assertEqual(summary["peak_width_q25"], 7.5)
        self.assertEqual(summary["peak_width_q75"], 12.5)

    def test_frip_counts_a_bedpe_pair_once_even_when_both_mates_overlap(self):
        result = calculate_frip([self.fragments[0]], [self.peaks[0]])

        self.assertEqual(result, {
            "total_fragments": 1,
            "fragments_in_peaks": 1,
            "frip": 1.0,
        })

    def test_frip_counts_fragment_once_when_it_overlaps_multiple_peaks(self):
        result = calculate_frip([self.fragments[1]], self.peaks[:2])

        self.assertEqual(result, {
            "total_fragments": 1,
            "fragments_in_peaks": 1,
            "frip": 1.0,
        })

    def test_frip_deduplicates_repeated_bedpe_records_by_fragment_name(self):
        result = calculate_frip([self.fragments[0], self.fragments[0]], [self.peaks[0]])

        self.assertEqual(result, {
            "total_fragments": 1,
            "fragments_in_peaks": 1,
            "frip": 1.0,
        })

    def test_frip_rejects_bedpe_records_without_a_fragment_name(self):
        for fragment in (
            ("chr1", 5, 10, "chr1", 20, 25),
            ("chr1", 5, 10, "chr1", 20, 25, ""),
            ("chr1", 5, 10, "chr1", 20, 25, "."),
        ):
            with self.subTest(fragment=fragment):
                with self.assertRaisesRegex(PeakQcError, "fragment name"):
                    calculate_frip([fragment], [self.peaks[0]])

    def test_named_bedpe_identity_uses_name_not_coordinates(self):
        result = calculate_frip(
            [
                ("chr1", 5, 10, "chr1", 20, 25, "duplicate_name"),
                ("chr1", 40, 45, "chr1", 50, 55, "duplicate_name"),
                ("chr1", 5, 10, "chr1", 20, 25, "distinct_name"),
            ],
            [self.peaks[0]],
        )

        self.assertEqual(result, {
            "total_fragments": 2,
            "fragments_in_peaks": 2,
            "frip": 1.0,
        })

    def test_overlap_analysis_avoids_all_fragment_peak_pairs(self):
        peaks = [Peak("chr1", index * 100, index * 100 + 10, f"peak_{index}") for index in range(100)]
        fragments = [
            Fragment(
                "chr1", index * 100, index * 100 + 2,
                "chr1", index * 100 + 8, index * 100 + 10,
                f"fragment_{index}",
            )
            for index in range(100)
        ]
        original_overlaps = peak_qc._overlaps

        with patch.object(peak_qc, "_overlaps", wraps=original_overlaps) as overlaps:
            rendered = peak_qc._render_outputs("sample", peaks, fragments)

        self.assertEqual(json.loads(rendered["json"])["fragments_in_peaks"], 100)
        self.assertLess(
            overlaps.call_count,
            300,
            "non-overlapping fragments and peaks must not trigger an all-pairs scan",
        )

    def test_write_qc_outputs_processes_fragments_incrementally(self):
        processed_names = []

        def guarded_fragments():
            yield ("chr1", 10, 12, "chr1", 18, 20, "fragment_one")
            if processed_names != ["fragment_one"]:
                raise AssertionError(
                    "fragment stream was consumed before overlap processing"
                )
            yield ("chr1", 30, 32, "chr1", 38, 40, "fragment_two")

        original_fragment_intervals = peak_qc._fragment_intervals

        def record_processed_fragment(fragment):
            processed_names.append(fragment.name)
            return original_fragment_intervals(fragment)

        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "streamed"
            with patch.object(
                peak_qc,
                "_fragment_intervals",
                side_effect=record_processed_fragment,
            ):
                write_qc_outputs(
                    prefix,
                    "sample",
                    [("chr1", 0, 50, "peak_a")],
                    guarded_fragments(),
                )

            payload = json.loads(Path(f"{prefix}.json").read_text(encoding="utf-8"))

        self.assertEqual(processed_names, ["fragment_one", "fragment_two"])
        self.assertEqual(payload["total_fragments"], 2)
        self.assertEqual(payload["fragments_in_peaks"], 2)

    def test_broadpeak_score_and_signal_value_are_summarized(self):
        peaks = read_peaks(ROOT / "tests" / "data" / "qc" / "scored.broadPeak")
        summary = summarize_peaks(peaks)

        self.assertEqual(peaks[0].score, 100)
        self.assertEqual(peaks[0].signal_value, 2.5)
        self.assertEqual(summary["peak_score_count"], 2)
        self.assertEqual(summary["peak_score_mean"], 200)
        self.assertEqual(summary["peak_score_q25"], 150)
        self.assertEqual(summary["peak_score_q75"], 250)
        self.assertEqual(summary["signal_value_count"], 2)
        self.assertEqual(summary["signal_value_mean"], 4.5)
        self.assertEqual(summary["signal_value_q25"], 3.5)
        self.assertEqual(summary["signal_value_q75"], 5.5)

    def test_broadpeak_score_and_signal_value_must_be_finite(self):
        for peak in (
            ("chr1", 10, 20, "bad_score", "NaN"),
            ("chr1", 10, 20, "bad_signal", "100", ".", "inf"),
            ("chr1", 10, 20, "bad_negative_signal", "100", ".", "-inf"),
            Peak("chr1", 10, 20, "constructed_nan", float("nan")),
            Peak("chr1", 10, 20, "constructed_inf", 100, float("inf")),
        ):
            with self.subTest(peak=peak):
                with self.assertRaisesRegex(PeakQcError, "finite"):
                    summarize_peaks([peak])

    def test_empty_peaks_have_zero_coverage_and_zero_frip(self):
        summary = summarize_peaks([])
        result = calculate_frip(self.fragments, [])

        self.assertEqual(summary["peak_count"], 0)
        self.assertEqual(summary["total_covered_bases"], 0)
        self.assertIsNone(summary["peak_width_min"])
        self.assertIsNone(summary["peak_width_q25"])
        self.assertEqual(result, {
            "total_fragments": 4,
            "fragments_in_peaks": 0,
            "frip": 0.0,
        })

    def test_rejects_intervals_with_non_positive_width(self):
        with self.assertRaisesRegex(PeakQcError, "positive width"):
            summarize_peaks([("chr1", 20, 20, "invalid")])


class PeakQcCliTests(unittest.TestCase):
    def make_workspace(self):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        return Path(temporary_directory.name)

    def run_cli(self, workspace, *arguments):
        return subprocess.run(
            [sys.executable, str(ROOT / "bin" / "peak_qc.py"), *arguments],
            cwd=workspace,
            capture_output=True,
            text=True,
            check=False,
        )

    def write_inputs(self, workspace, *, peaks=None, fragments=None):
        peaks_path = workspace / "sample.broadPeak"
        fragments_path = workspace / "sample.bedpe"
        peaks_path.write_text(
            "" if peaks is None else peaks,
            encoding="utf-8",
        )
        fragments_path.write_text(
            "" if fragments is None else fragments,
            encoding="utf-8",
        )
        return peaks_path, fragments_path

    def test_cli_writes_all_qc_outputs_and_per_peak_fragment_counts(self):
        workspace = self.make_workspace()
        peaks, fragments = self.write_inputs(
            workspace,
            peaks=(
                "chr1\t10\t20\tpeak_a\t100\t.\t2\t3\t1\n"
                "chr1\t15\t30\tpeak_b\t90\t.\t1\t2\t0\n"
            ),
            fragments=(
                "chr1\t5\t10\tchr1\t20\t25\tfragment_one\n"
                "chr1\t19\t21\tchr1\t21\t22\tfragment_two\n"
                "chr1\t30\t31\tchr1\t31\t32\tfragment_three\n"
            ),
        )
        prefix = workspace / "results" / "sample"

        completed = self.run_cli(
            workspace,
            "--fragments", str(fragments),
            "--peaks", str(peaks),
            "--sample-id", "sample",
            "--output-prefix", str(prefix),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(prefix.with_suffix(".json").read_text(encoding="utf-8"))
        self.assertEqual(payload["sample_id"], "sample")
        self.assertEqual(payload["total_fragments"], 3)
        self.assertEqual(payload["fragments_in_peaks"], 2)
        self.assertEqual(payload["frip"], 2 / 3)
        self.assertTrue(prefix.with_suffix(".tsv").is_file())
        histogram = prefix.with_name(prefix.name + ".width_histogram.tsv")
        self.assertEqual(
            histogram.read_text(encoding="utf-8").splitlines(),
            ["width\tpeak_count", "10\t1", "15\t1"],
        )
        per_peak = prefix.with_name(prefix.name + ".fragments_per_peak.tsv")
        self.assertEqual(
            per_peak.read_text(encoding="utf-8").splitlines(),
            [
                "chrom\tstart\tend\tpeak_name\twidth\tscore\tsignal_value\tfragment_count",
                "chr1\t10\t20\tpeak_a\t10\t100.0\t2.0\t2",
                "chr1\t15\t30\tpeak_b\t15\t90.0\t1.0\t2",
            ],
        )

    def test_cli_appends_extensions_to_a_dotted_output_prefix(self):
        workspace = self.make_workspace()
        peaks, fragments = self.write_inputs(
            workspace,
            peaks="chr1\t10\t20\tpeak_a\n",
            fragments="chr1\t5\t10\tchr1\t20\t25\tfragment_one\n",
        )
        prefix = workspace / "results" / "sample.qc"

        completed = self.run_cli(
            workspace,
            "--fragments", str(fragments), "--peaks", str(peaks),
            "--sample-id", "sample", "--output-prefix", str(prefix),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(Path(f"{prefix}.json").is_file())
        self.assertTrue(Path(f"{prefix}.tsv").is_file())
        self.assertTrue(Path(f"{prefix}.width_histogram.tsv").is_file())
        self.assertTrue(Path(f"{prefix}.fragments_per_peak.tsv").is_file())
        self.assertFalse(prefix.with_suffix(".json").exists())

    def test_cli_emits_broadpeak_score_and_signal_distributions(self):
        workspace = self.make_workspace()
        fragments = workspace / "sample.bedpe"
        fragments.write_text(
            "chr1\t5\t10\tchr1\t20\t25\tfragment_one\n",
            encoding="utf-8",
        )
        prefix = workspace / "scored"

        completed = self.run_cli(
            workspace,
            "--fragments", str(fragments),
            "--peaks", str(ROOT / "tests" / "data" / "qc" / "scored.broadPeak"),
            "--sample-id", "scored", "--output-prefix", str(prefix),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(Path(f"{prefix}.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["peak_score_mean"], 200)
        self.assertEqual(payload["signal_value_mean"], 4.5)
        self.assertIn("peak_score_q25\t150.0", Path(f"{prefix}.tsv").read_text(encoding="utf-8"))
        self.assertEqual(
            Path(f"{prefix}.fragments_per_peak.tsv").read_text(encoding="utf-8").splitlines(),
            [
                "chrom\tstart\tend\tpeak_name\twidth\tscore\tsignal_value\tfragment_count",
                "chr1\t10\t20\tpeak_a\t10\t100.0\t2.5\t1",
                "chr1\t30\t50\tpeak_b\t20\t300.0\t6.5\t0",
            ],
        )

    def test_cli_handles_empty_peak_file(self):
        workspace = self.make_workspace()
        peaks, fragments = self.write_inputs(
            workspace,
            peaks="",
            fragments="chr1\t0\t1\tchr1\t1\t2\tfragment_one\n",
        )
        prefix = workspace / "empty"

        completed = self.run_cli(
            workspace,
            "--fragments", str(fragments), "--peaks", str(peaks),
            "--sample-id", "empty", "--output-prefix", str(prefix),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(prefix.with_suffix(".json").read_text(encoding="utf-8"))
        self.assertEqual(payload["peak_count"], 0)
        self.assertEqual(payload["frip"], 0.0)

    def test_cli_rejects_malformed_input_without_replacing_existing_outputs(self):
        workspace = self.make_workspace()
        peaks, fragments = self.write_inputs(
            workspace,
            peaks="chr1\t20\t20\tinvalid\n",
            fragments="chr1\t0\t1\tchr1\t1\t2\tfragment_one\n",
        )
        prefix = workspace / "results" / "sample"
        existing = {}
        for path in (
            prefix.with_suffix(".json"),
            prefix.with_suffix(".tsv"),
            prefix.with_name(prefix.name + ".width_histogram.tsv"),
            prefix.with_name(prefix.name + ".fragments_per_peak.tsv"),
        ):
            path.parent.mkdir(parents=True, exist_ok=True)
            content = f"previous:{path.name}\n".encode()
            path.write_bytes(content)
            existing[path] = content

        completed = self.run_cli(
            workspace,
            "--fragments", str(fragments), "--peaks", str(peaks),
            "--sample-id", "sample", "--output-prefix", str(prefix),
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("positive width", completed.stderr)
        self.assertEqual({path: path.read_bytes() for path in existing}, existing)
        self.assertEqual(list(prefix.parent.glob(".peak-qc-*")), [])

    def test_cli_rejects_non_finite_broadpeak_values_before_output_publication(self):
        workspace = self.make_workspace()
        peaks, fragments = self.write_inputs(
            workspace,
            peaks="chr1\t10\t20\tinvalid\tNaN\t.\t2.5\t-1\t-1\n",
            fragments="chr1\t5\t10\tchr1\t20\t25\tfragment_one\n",
        )
        prefix = workspace / "results" / "sample"
        output = Path(f"{prefix}.json")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("previous JSON\n", encoding="utf-8")

        completed = self.run_cli(
            workspace,
            "--fragments", str(fragments), "--peaks", str(peaks),
            "--sample-id", "sample", "--output-prefix", str(prefix),
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("finite", completed.stderr)
        self.assertEqual(output.read_text(encoding="utf-8"), "previous JSON\n")
        self.assertFalse(Path(f"{prefix}.tsv").exists())

    def test_output_publication_rolls_back_if_one_replacement_fails(self):
        workspace = self.make_workspace()
        prefix = workspace / "results" / "sample"
        existing = {}
        for path in (
            prefix.with_suffix(".json"),
            prefix.with_suffix(".tsv"),
            prefix.with_name(prefix.name + ".width_histogram.tsv"),
            prefix.with_name(prefix.name + ".fragments_per_peak.tsv"),
        ):
            path.parent.mkdir(parents=True, exist_ok=True)
            content = f"previous:{path.name}\n".encode()
            path.write_bytes(content)
            existing[path] = content

        real_replace = peak_qc.os.replace
        calls = 0

        def fail_on_second_publication(source, destination):
            nonlocal calls
            calls += 1
            if calls == 4:
                raise OSError("simulated output replacement failure")
            return real_replace(source, destination)

        with patch.object(peak_qc.os, "replace", side_effect=fail_on_second_publication):
            with self.assertRaisesRegex(OSError, "simulated output replacement failure"):
                write_qc_outputs(
                    prefix,
                    "sample",
                    [("chr1", 10, 20, "peak_a")],
                    [("chr1", 5, 10, "chr1", 20, 25, "fragment_one")],
                )

        self.assertEqual({path: path.read_bytes() for path in existing}, existing)
        self.assertEqual(list(prefix.parent.glob(".peak-qc-*")), [])


if __name__ == "__main__":
    unittest.main()
