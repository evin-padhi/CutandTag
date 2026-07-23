import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

from manifest import ManifestValidationError, load_and_validate


FIELDNAMES = [
    "sample_id",
    "library_id",
    "input_group",
    "barcode",
    "assay_target",
    "is_control",
    "control_id",
    "expected_motif",
    "r1",
    "r2",
    "i2",
]


def valid_rows():
    rows = []
    for library_id, input_group, target_a, target_b, control_id in [
        ("NX702", "25K", "IgG", "CTCF", "NX702_IgG"),
        ("NX703", "50K", "IgG", "CTCF", "NX703_IgG"),
        ("NX701", "100K", "IgG", "CTCF", "NX701_IgG"),
        ("NX704", "25K", "GATA1", "RUNX1", "NX702_IgG"),
        ("NX705", "50K", "GATA1", "RUNX1", "NX703_IgG"),
        ("NX706", "100K", "GATA1", "RUNX1", "NX701_IgG"),
    ]:
        for barcode, target in [("TATAGCCT", target_a), ("ATAGAGGC", target_b)]:
            control = target == "IgG"
            rows.append(
                {
                    "sample_id": f"{library_id}_{target}",
                    "library_id": library_id,
                    "input_group": input_group,
                    "barcode": barcode,
                    "assay_target": target,
                    "is_control": str(control).lower(),
                    "control_id": "" if control else control_id,
                    "expected_motif": "" if control else target,
                    "r1": f"data/{library_id}_R1.fastq.gz",
                    "r2": f"data/{library_id}_R2.fastq.gz",
                    "i2": f"data/{library_id}_I2.fastq.gz",
                }
            )
    return rows


def write_manifest(tmp_path, rows):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    for library_id in {row["library_id"] for row in rows}:
        for read in ("R1", "R2", "I2"):
            (data_dir / f"{library_id}_{read}.fastq.gz").touch()

    manifest = tmp_path / "samples.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return manifest


def write_manifest_without_fastqs(tmp_path, rows):
    manifest = tmp_path / "samples.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return manifest


class ManifestTests(unittest.TestCase):
    def temporary_manifest(self, rows):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        temp_path = Path(temporary_directory.name)
        return temp_path, write_manifest(temp_path, rows)

    def test_valid_six_library_manifest_normalizes_records_and_relative_paths(self):
        temp_path, manifest = self.temporary_manifest(valid_rows())

        records = load_and_validate(manifest)

        self.assertEqual(len(records), 12)
        self.assertIs(records[0]["is_control"], True)
        self.assertIs(records[1]["is_control"], False)
        self.assertIsNone(records[0]["control_id"])
        self.assertIsNone(records[0]["expected_motif"])
        self.assertEqual(records[1]["r1"], str((temp_path / "data/NX702_R1.fastq.gz").resolve()))
        self.assertEqual(
            {record["library_id"] for record in records},
            {"NX701", "NX702", "NX703", "NX704", "NX705", "NX706"},
        )

    def test_duplicate_sample_ids_are_rejected(self):
        rows = valid_rows()
        rows[1]["sample_id"] = rows[0]["sample_id"]

        _, manifest = self.temporary_manifest(rows)
        with self.assertRaisesRegex(ManifestValidationError, "duplicate sample_id"):
            load_and_validate(manifest)

    def test_sample_ids_must_be_safe_filename_components(self):
        invalid_identifiers = (
            "/absolute",
            "nested/sample",
            r"nested\sample",
            "has whitespace",
            ".leading-dot",
            "-leading-dash",
            "shell;command",
            "shell$(command)",
        )
        for invalid_identifier in invalid_identifiers:
            with self.subTest(identifier=invalid_identifier):
                rows = valid_rows()
                rows[1]["sample_id"] = invalid_identifier
                with tempfile.TemporaryDirectory() as temporary_directory:
                    manifest = write_manifest_without_fastqs(
                        Path(temporary_directory), rows
                    )
                    with self.assertRaisesRegex(
                        ManifestValidationError, "invalid sample_id"
                    ):
                        load_and_validate(manifest, check_files=False)

    def test_library_ids_must_be_safe_filename_components(self):
        invalid_identifiers = (
            "/absolute",
            "nested/library",
            r"nested\library",
            "has whitespace",
            ".leading-dot",
            "-leading-dash",
            "shell;command",
            "shell${command}",
        )
        for invalid_identifier in invalid_identifiers:
            with self.subTest(identifier=invalid_identifier):
                rows = valid_rows()
                rows[0]["library_id"] = invalid_identifier
                with tempfile.TemporaryDirectory() as temporary_directory:
                    manifest = write_manifest_without_fastqs(
                        Path(temporary_directory), rows
                    )
                    with self.assertRaisesRegex(
                        ManifestValidationError, "invalid library_id"
                    ):
                        load_and_validate(manifest, check_files=False)

    def test_safe_identifier_punctuation_is_accepted(self):
        rows = valid_rows()
        rows[1]["sample_id"] = "A.sample-1_target"
        rows[0]["library_id"] = "A.library-1_control"
        rows[1]["library_id"] = "A.library-1_control"

        with tempfile.TemporaryDirectory() as temporary_directory:
            manifest = write_manifest_without_fastqs(
                Path(temporary_directory), rows
            )
            records = load_and_validate(manifest, check_files=False)

        self.assertEqual(records[1]["sample_id"], "A.sample-1_target")
        self.assertEqual(records[0]["library_id"], "A.library-1_control")

    def test_inconsistent_fastq_paths_within_a_library_are_rejected(self):
        rows = valid_rows()
        rows[1]["r1"] = "data/other_R1.fastq.gz"
        temp_path, manifest = self.temporary_manifest(rows)
        (temp_path / "data/other_R1.fastq.gz").touch()

        with self.assertRaisesRegex(ManifestValidationError, "inconsistent FASTQ paths"):
            load_and_validate(manifest)

    def test_unequal_barcode_lengths_are_rejected(self):
        rows = valid_rows()
        rows[0]["barcode"] = "TATAGCC"

        _, manifest = self.temporary_manifest(rows)
        with self.assertRaisesRegex(ManifestValidationError, "equal length"):
            load_and_validate(manifest)

    def test_duplicate_barcodes_within_a_library_are_rejected(self):
        rows = valid_rows()
        rows[1]["barcode"] = rows[0]["barcode"]

        _, manifest = self.temporary_manifest(rows)
        with self.assertRaisesRegex(ManifestValidationError, "duplicate barcode"):
            load_and_validate(manifest)

    def test_missing_controls_are_rejected(self):
        rows = valid_rows()
        rows[1]["control_id"] = "UNKNOWN_IgG"

        _, manifest = self.temporary_manifest(rows)
        with self.assertRaisesRegex(ManifestValidationError, "missing control"):
            load_and_validate(manifest)

    def test_non_igg_controls_are_rejected(self):
        rows = valid_rows()
        rows[0]["assay_target"] = "CTCF"

        _, manifest = self.temporary_manifest(rows)
        with self.assertRaisesRegex(ManifestValidationError, "must have assay_target IgG"):
            load_and_validate(manifest)

    def test_targets_and_controls_must_share_input_group(self):
        rows = valid_rows()
        rows[6]["input_group"] = "50K"

        _, manifest = self.temporary_manifest(rows)
        with self.assertRaisesRegex(ManifestValidationError, "input_group"):
            load_and_validate(manifest)

    def test_cli_writes_normalized_json(self):
        temp_path, manifest = self.temporary_manifest(valid_rows())
        output = temp_path / "normalized.json"

        completed = subprocess.run(
            [
                sys.executable,
                str(ROOT / "bin/manifest.py"),
                "validate",
                "--input",
                str(manifest),
                "--output",
                str(output),
            ],
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(json.loads(output.read_text())[0]["sample_id"], "NX702_IgG")

    def test_cli_keeps_file_checks_strict_by_default(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "normalized.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "bin/manifest.py"),
                    "validate",
                    "--input",
                    str(ROOT / "assets/samples.example.csv"),
                    "--output",
                    str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
            )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("file does not exist", completed.stderr)

    def test_cli_skip_file_checks_validates_the_committed_example(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "normalized.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "bin/manifest.py"),
                    "validate",
                    "--input",
                    str(ROOT / "assets/samples.example.csv"),
                    "--output",
                    str(output),
                    "--skip-file-checks",
                ],
                check=False,
                capture_output=True,
                text=True,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(len(json.loads(output.read_text())), 12)


if __name__ == "__main__":
    unittest.main()
