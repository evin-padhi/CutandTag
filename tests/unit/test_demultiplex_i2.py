import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

from demultiplex_i2 import assign_barcode, hamming


def fastq_record(identifier, sequence):
    return f"@{identifier}\n{sequence}\n+\n{'I' * len(sequence)}\n"


def write_fastq(path, records, *, compressed=False):
    opener = gzip.open if compressed else open
    mode = "wt"
    with opener(path, mode, encoding="utf-8") as handle:
        for identifier, sequence in records:
            handle.write(fastq_record(identifier, sequence))


def read_fastq_headers(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        lines = handle.readlines()
    return [line.rstrip() for line in lines[::4]]


def write_previous_outputs(outdir, metrics, sample_ids):
    outdir.mkdir(parents=True, exist_ok=True)
    previous = {}
    for sample_id in sample_ids:
        for read in ("R1", "R2"):
            path = outdir / f"{sample_id}_{read}.fastq.gz"
            content = f"previous {sample_id} {read}".encode()
            path.write_bytes(content)
            previous[path] = content
    metrics_tsv = metrics.with_suffix(".tsv")
    for path, content in (
        (metrics, b'{"previous": true}\n'),
        (metrics_tsv, b"metric\tvalue\nprevious\t1\n"),
    ):
        path.write_bytes(content)
        previous[path] = content
    return previous


class DemultiplexI2UnitTests(unittest.TestCase):
    expected = {"sample_A": "AAAAAAAA", "sample_B": "CCCCCCCC"}

    def test_hamming_and_exact_assignments(self):
        self.assertEqual(hamming("AAAA", "AAAT"), 1)
        self.assertEqual(assign_barcode("AAAAAAAA", self.expected, 0), ("assigned", "sample_A"))
        self.assertEqual(assign_barcode("CCCCCCCC", self.expected, 0), ("assigned", "sample_B"))

    def test_mismatch_override_selects_the_unique_closest_barcode(self):
        expected = {"sample_A": "AAAAAAAA", "sample_B": "AAAACCCC"}

        self.assertEqual(assign_barcode("AAAAAAAT", expected, 0), ("unassigned", None))
        self.assertEqual(assign_barcode("AAAAAAAT", expected, 1), ("assigned", "sample_A"))
        self.assertEqual(assign_barcode("AAAACCCA", expected, 1), ("assigned", "sample_B"))

    def test_equal_closest_barcodes_are_ambiguous(self):
        expected = {"sample_A": "AAAAAAAA", "sample_B": "AAAAAAAT"}

        self.assertEqual(assign_barcode("AAAAAAAG", expected, 1), ("ambiguous", None))

    def test_reads_outside_the_mismatch_threshold_are_unassigned(self):
        self.assertEqual(assign_barcode("GGGGGGGG", self.expected, 2), ("unassigned", None))


class DemultiplexI2CliTests(unittest.TestCase):
    def make_workspace(self):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        return Path(temporary_directory.name)

    def run_cli(self, workspace, *extra_arguments):
        return subprocess.run(
            [
                sys.executable,
                str(ROOT / "bin/demultiplex_i2.py"),
                *extra_arguments,
            ],
            capture_output=True,
            text=True,
            check=False,
            cwd=workspace,
        )

    def test_cli_splits_exact_assignments_and_writes_metrics(self):
        workspace = self.make_workspace()
        r1 = workspace / "reads_R1.fastq"
        r2 = workspace / "reads_R2.fastq"
        i2 = workspace / "reads_I2.fastq"
        write_fastq(r1, [("read-a/1 1:N:0:1", "ACGT"), ("read-b/1", "TGCA"), ("read-c/1", "GCGC")])
        write_fastq(r2, [("read-a/2 2:N:0:1", "TTTT"), ("read-b/2", "CCCC"), ("read-c/2", "AAAA")])
        write_fastq(i2, [("read-a 3:N:0:1", "AAAAAAAA"), ("read-b", "CCCCCCCC"), ("read-c", "GGGGGGGG")])
        outdir = workspace / "out"
        metrics = workspace / "metrics.json"

        completed = self.run_cli(
            workspace,
            "--r1", str(r1), "--r2", str(r2), "--i2", str(i2),
            "--sample", "sample_A=AAAAAAAA", "--sample", "sample_B=CCCCCCCC",
            "--outdir", str(outdir), "--metrics", str(metrics),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(read_fastq_headers(outdir / "sample_A_R1.fastq.gz"), ["@read-a/1 1:N:0:1"])
        self.assertEqual(read_fastq_headers(outdir / "sample_B_R2.fastq.gz"), ["@read-b/2"])
        payload = json.loads(metrics.read_text())
        self.assertEqual(payload["total_reads"], 3)
        self.assertEqual(payload["assigned_reads"], 2)
        self.assertEqual(payload["unassigned_reads"], 1)
        self.assertEqual(payload["ambiguous_reads"], 0)
        self.assertEqual(payload["assignment_counts"], {"sample_A": 1, "sample_B": 1})
        self.assertEqual(payload["observed_i2_counts"], {"AAAAAAAA": 1, "CCCCCCCC": 1, "GGGGGGGG": 1})
        self.assertEqual(payload["assigned_fraction"], 2 / 3)
        self.assertTrue((workspace / "metrics.tsv").is_file())

    def test_cli_accepts_a_mixture_of_gzip_and_plain_fastqs(self):
        workspace = self.make_workspace()
        r1 = workspace / "reads_R1.fastq.gz"
        r2 = workspace / "reads_R2.fastq"
        i2 = workspace / "reads_I2.fastq.gz"
        write_fastq(r1, [("read-a/1", "ACGT")], compressed=True)
        write_fastq(r2, [("read-a/2", "TTTT")])
        write_fastq(i2, [("read-a", "AAAAAAAA")], compressed=True)

        completed = self.run_cli(
            workspace,
            "--r1", str(r1), "--r2", str(r2), "--i2", str(i2),
            "--sample", "sample_A=AAAAAAAA", "--outdir", str(workspace / "out"),
            "--metrics", str(workspace / "metrics.json"),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(read_fastq_headers(workspace / "out/sample_A_R1.fastq.gz"), ["@read-a/1"])

    def test_cli_accepts_synchronized_read_number_suffixes(self):
        workspace = self.make_workspace()
        r1 = workspace / "reads_R1.fastq"
        r2 = workspace / "reads_R2.fastq"
        i2 = workspace / "reads_I2.fastq"
        write_fastq(r1, [("read-a/1", "ACGT")])
        write_fastq(r2, [("read-a/2", "TTTT")])
        write_fastq(i2, [("read-a/3", "AAAAAAAA")])

        completed = self.run_cli(
            workspace,
            "--r1", str(r1), "--r2", str(r2), "--i2", str(i2),
            "--sample", "sample_A=AAAAAAAA", "--outdir", str(workspace / "out"),
            "--metrics", str(workspace / "metrics.json"),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(read_fastq_headers(workspace / "out/sample_A_R1.fastq.gz"), ["@read-a/1"])

    def test_cli_rejects_truncated_fastq_records(self):
        workspace = self.make_workspace()
        r1 = workspace / "reads_R1.fastq"
        r2 = workspace / "reads_R2.fastq"
        i2 = workspace / "reads_I2.fastq"
        r1.write_text("@read-a/1\nACGT\n+\n")
        write_fastq(r2, [("read-a/2", "TTTT")])
        write_fastq(i2, [("read-a", "AAAAAAAA")])
        outdir = workspace / "out"
        metrics = workspace / "metrics.json"
        previous = write_previous_outputs(outdir, metrics, ["sample_A"])

        completed = self.run_cli(
            workspace,
            "--r1", str(r1), "--r2", str(r2), "--i2", str(i2),
            "--sample", "sample_A=AAAAAAAA", "--outdir", str(outdir),
            "--metrics", str(metrics),
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("truncated FASTQ", completed.stderr)
        self.assertEqual({path: path.read_bytes() for path in previous}, previous)
        self.assertEqual(list(outdir.glob(".demultiplex-*")), [])

    def test_cli_rejects_normalized_read_identifier_mismatches(self):
        workspace = self.make_workspace()
        r1 = workspace / "reads_R1.fastq"
        r2 = workspace / "reads_R2.fastq"
        i2 = workspace / "reads_I2.fastq"
        write_fastq(r1, [("read-a/1", "ACGT")])
        write_fastq(r2, [("read-other/2", "TTTT")])
        write_fastq(i2, [("read-a", "AAAAAAAA")])

        completed = self.run_cli(
            workspace,
            "--r1", str(r1), "--r2", str(r2), "--i2", str(i2),
            "--sample", "sample_A=AAAAAAAA", "--outdir", str(workspace / "out"),
            "--metrics", str(workspace / "metrics.json"),
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("read identifiers differ", completed.stderr)

    def test_cli_rejects_unequal_fastq_record_counts(self):
        workspace = self.make_workspace()
        r1 = workspace / "reads_R1.fastq"
        r2 = workspace / "reads_R2.fastq"
        i2 = workspace / "reads_I2.fastq"
        write_fastq(r1, [("read-a/1", "ACGT"), ("read-b/1", "TGCA")])
        write_fastq(r2, [("read-a/2", "TTTT")])
        write_fastq(i2, [("read-a", "AAAAAAAA")])

        completed = self.run_cli(
            workspace,
            "--r1", str(r1), "--r2", str(r2), "--i2", str(i2),
            "--sample", "sample_A=AAAAAAAA", "--outdir", str(workspace / "out"),
            "--metrics", str(workspace / "metrics.json"),
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("record counts differ", completed.stderr)

    def test_empty_derived_sample_fails_unless_explicitly_allowed(self):
        workspace = self.make_workspace()
        r1 = workspace / "reads_R1.fastq"
        r2 = workspace / "reads_R2.fastq"
        i2 = workspace / "reads_I2.fastq"
        write_fastq(r1, [("read-a/1", "ACGT")])
        write_fastq(r2, [("read-a/2", "TTTT")])
        write_fastq(i2, [("read-a", "AAAAAAAA")])
        outdir = workspace / "out"
        metrics = workspace / "metrics.json"
        previous = write_previous_outputs(outdir, metrics, ["sample_A", "sample_B"])
        common_arguments = (
            "--r1", str(r1), "--r2", str(r2), "--i2", str(i2),
            "--sample", "sample_A=AAAAAAAA", "--sample", "sample_B=CCCCCCCC",
            "--outdir", str(outdir), "--metrics", str(metrics),
        )

        rejected = self.run_cli(workspace, *common_arguments)

        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn("empty derived samples", rejected.stderr)
        self.assertEqual({path: path.read_bytes() for path in previous}, previous)
        self.assertEqual(list(outdir.glob(".demultiplex-*")), [])

        allowed = self.run_cli(workspace, *common_arguments, "--allow-empty")
        self.assertEqual(allowed.returncode, 0, allowed.stderr)

    def test_failed_run_preserves_every_preexisting_derived_output_and_metric(self):
        workspace = self.make_workspace()
        r1 = workspace / "reads_R1.fastq"
        r2 = workspace / "reads_R2.fastq"
        i2 = workspace / "reads_I2.fastq"
        write_fastq(r1, [("read-a/1", "ACGT")])
        write_fastq(r2, [("wrong-read/2", "TTTT")])
        write_fastq(i2, [("read-a", "AAAAAAAA")])
        outdir = workspace / "out"
        metrics = workspace / "metrics.json"
        previous = write_previous_outputs(outdir, metrics, ["sample_A", "sample_B"])

        completed = self.run_cli(
            workspace,
            "--r1", str(r1), "--r2", str(r2), "--i2", str(i2),
            "--sample", "sample_A=AAAAAAAA", "--sample", "sample_B=CCCCCCCC",
            "--outdir", str(outdir), "--metrics", str(metrics),
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("read identifiers differ", completed.stderr)
        self.assertEqual({path: path.read_bytes() for path in previous}, previous)
        self.assertEqual(list(outdir.glob(".demultiplex-*")), [])

    def test_cli_rejects_an_invalid_mismatch_limit(self):
        workspace = self.make_workspace()
        for read in ("R1", "R2", "I2"):
            write_fastq(workspace / f"reads_{read}.fastq", [(f"read-a/{1 if read == 'R1' else 2}", "AAAAAAAA")])

        completed = self.run_cli(
            workspace,
            "--r1", str(workspace / "reads_R1.fastq"), "--r2", str(workspace / "reads_R2.fastq"),
            "--i2", str(workspace / "reads_I2.fastq"), "--sample", "sample_A=AAAAAAAA",
            "--max-mismatches", "8", "--outdir", str(workspace / "out"),
            "--metrics", str(workspace / "metrics.json"),
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("smaller than barcode length", completed.stderr)

    def test_cli_rejects_a_negative_mismatch_limit(self):
        workspace = self.make_workspace()
        for read in ("R1", "R2", "I2"):
            write_fastq(workspace / f"reads_{read}.fastq", [("read-a", "AAAAAAAA")])

        completed = self.run_cli(
            workspace,
            "--r1", str(workspace / "reads_R1.fastq"), "--r2", str(workspace / "reads_R2.fastq"),
            "--i2", str(workspace / "reads_I2.fastq"), "--sample", "sample_A=AAAAAAAA",
            "--max-mismatches", "-1", "--outdir", str(workspace / "out"),
            "--metrics", str(workspace / "metrics.json"),
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("non-negative", completed.stderr)


if __name__ == "__main__":
    unittest.main()
