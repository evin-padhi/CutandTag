import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

import motif_qc
from motif_qc import MotifQcError, read_ame, read_fimo, summarize_expected_motif, write_motif_outputs


class MotifQcUnitTests(unittest.TestCase):
    def setUp(self):
        self.ame = read_ame(ROOT / "tests" / "data" / "motifs" / "ame.tsv")
        self.fimo = read_fimo(ROOT / "tests" / "data" / "motifs" / "fimo.tsv")

    def test_regex_selects_expected_motifs_and_ranks_by_adjusted_p_value(self):
        summary, positions = summarize_expected_motif(
            self.ame, self.fimo, r"^CTCF", window=200
        )

        self.assertEqual(summary["status"], "pass")
        self.assertEqual(summary["matching_motif_count"], 2)
        self.assertEqual(summary["best_motif_id"], "MA0139.2")
        self.assertEqual(summary["best_motif_rank"], 2)
        self.assertEqual(summary["best_adjusted_p_value"], 0.002)
        self.assertEqual(summary["best_effect"], 4.5)
        self.assertEqual(len(positions), 3)

    def test_multiple_expected_motifs_use_unique_peak_coverage(self):
        summary, _ = summarize_expected_motif(self.ame, self.fimo, r"^CTCF", window=200)

        self.assertEqual(summary["peak_count"], 4)
        self.assertEqual(summary["expected_hit_count"], 3)
        self.assertEqual(summary["peaks_with_expected_motif"], 2)
        self.assertEqual(summary["peak_fraction_with_expected_motif"], 0.5)

    def test_hit_positions_are_relative_to_the_window_center(self):
        _, positions = summarize_expected_motif(self.ame, self.fimo, r"^CTCF", window=200)

        self.assertEqual(
            [(row["sequence_name"], row["relative_center_position"]) for row in positions],
            [("peak_1", 0.0), ("peak_1", 0.0), ("peak_2", -95.0)],
        )
        self.assertEqual(sum(row["is_central"] for row in positions), 2)

    def test_no_matching_database_motif_is_reported_explicitly(self):
        summary, positions = summarize_expected_motif(self.ame, self.fimo, r"^GATA1$", window=200)

        self.assertEqual(summary["status"], "motif_not_found")
        self.assertEqual(summary["matching_motif_count"], 0)
        self.assertEqual(positions, [])

    def test_empty_ame_and_fimo_results_are_no_peaks(self):
        summary, positions = summarize_expected_motif([], [], r"CTCF", window=200)

        self.assertEqual(summary["status"], "no_peaks")
        self.assertEqual(positions, [])

    def test_header_only_ame_report_is_not_mistaken_for_no_peaks(self):
        workspace = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(workspace))
        ame = workspace / "ame.tsv"
        fimo = workspace / "fimo.tsv"
        ame.write_text(
            "rank\tmotif_ID\tmotif_Alt_ID\tadj_p-value\tpos\n",
            encoding="utf-8",
        )
        fimo.write_text(
            "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\n",
            encoding="utf-8",
        )

        summary, positions = summarize_expected_motif(
            read_ame(ame), read_fimo(fimo), r"^CTCF$", window=200
        )

        self.assertEqual(summary["status"], "not_significant")
        self.assertEqual(positions, [])

    def test_non_significant_expected_motif_has_distinct_status(self):
        summary, _ = summarize_expected_motif(self.ame, self.fimo, r"^SPI1$", window=200)

        self.assertEqual(summary["status"], "not_significant")
        self.assertEqual(summary["best_adjusted_p_value"], 0.2)

    def test_unmatched_default_ame_report_is_conservatively_not_significant(self):
        workspace = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(workspace))
        ame = workspace / "ame.tsv"
        ame.write_text(
            "rank\tmotif_ID\tmotif_Alt_ID\tadj_p-value\tpos\n"
            "1\tMA0476.1\tFOS\t0.001\t4\n",
            encoding="utf-8",
        )

        summary, _ = summarize_expected_motif(read_ame(ame), [], r"^CTCF$", window=200)

        self.assertEqual(summary["status"], "not_significant")

    def test_modern_ame_e_value_is_supported_as_the_adjusted_statistic(self):
        workspace = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(workspace))
        ame = workspace / "ame.tsv"
        ame.write_text(
            "# motif_qc_complete_database=true\n"
            "rank\tmotif_ID\tmotif_Alt_ID\tp-value\tE-value\tpos\tneg\n"
            "1\tMA0139.1\tCTCF\t0.0001\t2.0\t4\t4\n",
            encoding="utf-8",
        )

        records = read_ame(ame)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].adjusted_p_value, 2.0)
        self.assertTrue(records.complete_database)

    def test_modern_fisher_ame_derives_hit_rate_enrichment_effect(self):
        workspace = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(workspace))
        ame = workspace / "ame.tsv"
        ame.write_text(
            "rank\tmotif_ID\tmotif_Alt_ID\tp-value\tE-value\tpos\tneg\tTP\t%TP\tFP\t%FP\n"
            "1\tMA0139.1\tCTCF\t0.0001\t0.002\t20\t40\t10\t50\t5\t12.5\n",
            encoding="utf-8",
        )

        records = read_ame(ame)
        summary, _ = summarize_expected_motif(records, [], r"^CTCF$", window=200)

        self.assertEqual(records[0].effect, 4.0)
        self.assertEqual(summary["best_effect"], 4.0)

    def test_readers_reject_invalid_probability_statistics(self):
        workspace = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(workspace))
        for value, expression in (("NaN", "finite"), ("-0.1", "between 0 and 1"), ("1.1", "between 0 and 1")):
            with self.subTest(value=value):
                ame = workspace / f"bad_ame_{value}.tsv"
                ame.write_text(
                    f"motif_ID\tmotif_Alt_ID\tadj_p-value\tpos\nMA1\tCTCF\t{value}\t1\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(MotifQcError, expression):
                    read_ame(ame)
        fimo = workspace / "bad_fimo.tsv"
        fimo.write_text(
            "motif_id\tsequence_name\tstart\tstop\tp-value\tq-value\n"
            "MA1\tpeak\t1\t2\t0.5\t1.2\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(MotifQcError, "between 0 and 1"):
            read_fimo(fimo)


class MotifQcCliTests(unittest.TestCase):
    def make_workspace(self):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        return Path(temporary_directory.name)

    def run_cli(self, workspace, *arguments):
        return subprocess.run(
            [sys.executable, str(ROOT / "bin" / "motif_qc.py"), *arguments],
            cwd=workspace,
            capture_output=True,
            text=True,
            check=False,
        )

    def test_cli_writes_json_tsv_and_centered_hit_positions(self):
        workspace = self.make_workspace()
        prefix = workspace / "results" / "sample.qc"

        completed = self.run_cli(
            workspace,
            "--ame", str(ROOT / "tests" / "data" / "motifs" / "ame.tsv"),
            "--fimo", str(ROOT / "tests" / "data" / "motifs" / "fimo.tsv"),
            "--expected-motif", r"^CTCF", "--window", "200",
            "--output-prefix", str(prefix),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(Path(f"{prefix}.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "pass")
        self.assertEqual(payload["peak_fraction_with_expected_motif"], 0.5)
        self.assertEqual(
            Path(f"{prefix}.motif_hit_positions.tsv").read_text(encoding="utf-8").splitlines(),
            [
                "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tp_value\tq_value\trelative_center_position\tis_central",
                "MA0139.1\tCTCF\tpeak_1\t96\t105\t+\t0.0001\t0.001\t0.0\ttrue",
                "MA0139.2\tCTCF\tpeak_1\t96\t105\t-\t0.0002\t0.002\t0.0\ttrue",
                "MA0139.2\tCTCF\tpeak_2\t1\t10\t+\t0.02\t0.1\t-95.0\tfalse",
            ],
        )
        self.assertTrue(Path(f"{prefix}.tsv").is_file())

    def test_invalid_input_preserves_existing_outputs(self):
        workspace = self.make_workspace()
        bad_ame = workspace / "bad.tsv"
        bad_ame.write_text("motif_ID\tadj_p-value\nMA1\tinf\n", encoding="utf-8")
        prefix = workspace / "results" / "sample"
        outputs = [
            Path(f"{prefix}.json"), Path(f"{prefix}.tsv"),
            Path(f"{prefix}.motif_hit_positions.tsv"),
        ]
        before = {}
        for output in outputs:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(f"old:{output.name}\n", encoding="utf-8")
            before[output] = output.read_bytes()

        completed = self.run_cli(
            workspace, "--ame", str(bad_ame),
            "--fimo", str(ROOT / "tests" / "data" / "motifs" / "fimo.tsv"),
            "--expected-motif", "CTCF", "--window", "200", "--output-prefix", str(prefix),
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("finite", completed.stderr)
        self.assertEqual({path: path.read_bytes() for path in outputs}, before)
        self.assertEqual(list(prefix.parent.glob(".motif-qc-*")), [])

    def test_publication_rolls_back_when_a_replacement_fails(self):
        workspace = self.make_workspace()
        prefix = workspace / "results" / "sample"
        outputs = [
            Path(f"{prefix}.json"), Path(f"{prefix}.tsv"),
            Path(f"{prefix}.motif_hit_positions.tsv"),
        ]
        before = {}
        for output in outputs:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(f"old:{output.name}\n", encoding="utf-8")
            before[output] = output.read_bytes()
        summary, positions = summarize_expected_motif(
            read_ame(ROOT / "tests" / "data" / "motifs" / "ame.tsv"),
            read_fimo(ROOT / "tests" / "data" / "motifs" / "fimo.tsv"),
            r"^CTCF", 200,
        )
        real_replace = motif_qc.os.replace
        calls = 0

        def fail_on_second_publish(source, destination):
            nonlocal calls
            calls += 1
            if calls == 4:
                raise OSError("simulated publication failure")
            return real_replace(source, destination)

        with patch.object(motif_qc.os, "replace", side_effect=fail_on_second_publish):
            with self.assertRaisesRegex(OSError, "simulated publication failure"):
                write_motif_outputs(prefix, summary, positions)

        self.assertEqual({path: path.read_bytes() for path in outputs}, before)
        self.assertEqual(list(prefix.parent.glob(".motif-qc-*")), [])


if __name__ == "__main__":
    unittest.main()
