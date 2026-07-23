#!/usr/bin/env python3
"""Validate and normalize bulk nano-CUT&Tag sample manifests."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import tempfile
from pathlib import Path


REQUIRED_COLUMNS = (
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
)
REQUIRED_VALUES = (
    "sample_id",
    "library_id",
    "input_group",
    "barcode",
    "assay_target",
    "is_control",
    "r1",
    "r2",
    "i2",
)
FASTQ_COLUMNS = ("r1", "r2", "i2")
ASSAY_TARGETS = {"IgG", "CTCF", "GATA1", "RUNX1"}
TRUE_VALUES = {"true", "1", "yes"}
FALSE_VALUES = {"false", "0", "no"}
IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class ManifestValidationError(ValueError):
    """Raised when a manifest violates the pipeline input contract."""


def _nonempty(value: object) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    return normalized or None


def _normalize_boolean(value: str, row_number: int) -> bool:
    normalized = value.lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    raise ManifestValidationError(
        f"row {row_number}: is_control must be a boolean, got {value!r}"
    )


def _resolve_fastq(
    value: str,
    manifest_directory: Path,
    row_number: int,
    column: str,
    check_files: bool,
) -> str:
    path = Path(value)
    if not path.is_absolute():
        path = manifest_directory / path
    path = path.resolve()
    if check_files and not path.is_file():
        raise ManifestValidationError(
            f"row {row_number}: {column} file does not exist: {path}"
        )
    return str(path)


def _read_rows(path: Path, check_files: bool) -> list[dict[str, object]]:
    if not path.is_file():
        raise ManifestValidationError(f"manifest file does not exist: {path}")

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        missing_columns = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
        if missing_columns:
            raise ManifestValidationError(
                f"manifest is missing required columns: {', '.join(missing_columns)}"
            )

        records: list[dict[str, object]] = []
        for row_number, row in enumerate(reader, start=2):
            values = {column: _nonempty(row.get(column)) for column in REQUIRED_COLUMNS}
            missing_values = [column for column in REQUIRED_VALUES if values[column] is None]
            if missing_values:
                raise ManifestValidationError(
                    f"row {row_number}: required values are blank: {', '.join(missing_values)}"
                )

            for column in ("sample_id", "library_id"):
                identifier = values[column]
                assert isinstance(identifier, str)
                if IDENTIFIER_PATTERN.fullmatch(identifier) is None:
                    raise ManifestValidationError(
                        f"row {row_number}: invalid {column} {identifier!r}; "
                        "use only letters, digits, dot, underscore, or dash, "
                        "and begin with a letter or digit"
                    )

            assay_target = values["assay_target"]
            assert isinstance(assay_target, str)
            if assay_target not in ASSAY_TARGETS:
                raise ManifestValidationError(
                    f"row {row_number}: unsupported assay_target {assay_target!r}"
                )

            is_control_value = values["is_control"]
            assert isinstance(is_control_value, str)
            is_control = _normalize_boolean(is_control_value, row_number)
            normalized: dict[str, object] = {
                "sample_id": values["sample_id"],
                "library_id": values["library_id"],
                "input_group": values["input_group"],
                "barcode": values["barcode"],
                "assay_target": assay_target,
                "is_control": is_control,
                "control_id": values["control_id"],
                "expected_motif": values["expected_motif"],
            }
            for column in FASTQ_COLUMNS:
                value = values[column]
                assert isinstance(value, str)
                normalized[column] = _resolve_fastq(
                    value, path.parent, row_number, column, check_files
                )
            normalized["_row_number"] = row_number
            records.append(normalized)

    if not records:
        raise ManifestValidationError("manifest contains no sample rows")
    return records


def _validate_relationships(records: list[dict[str, object]]) -> None:
    samples: dict[str, dict[str, object]] = {}
    library_paths: dict[str, tuple[str, str, str]] = {}
    library_barcodes: dict[str, dict[str, str]] = {}
    barcode_lengths = {len(str(record["barcode"])) for record in records}
    if len(barcode_lengths) != 1 or 0 in barcode_lengths:
        raise ManifestValidationError("all barcodes must be non-empty and equal length")

    for record in records:
        row_number = record["_row_number"]
        sample_id = str(record["sample_id"])
        if sample_id in samples:
            raise ManifestValidationError(f"row {row_number}: duplicate sample_id {sample_id!r}")
        samples[sample_id] = record

        library_id = str(record["library_id"])
        paths = tuple(str(record[column]) for column in FASTQ_COLUMNS)
        previous_paths = library_paths.setdefault(library_id, paths)
        if paths != previous_paths:
            raise ManifestValidationError(
                f"row {row_number}: library {library_id!r} has inconsistent FASTQ paths"
            )
        barcode = str(record["barcode"])
        previous_sample = library_barcodes.setdefault(library_id, {}).get(barcode)
        if previous_sample is not None:
            raise ManifestValidationError(
                f"row {row_number}: library {library_id!r} has duplicate barcode "
                f"{barcode!r} for samples {previous_sample!r} and {sample_id!r}"
            )
        library_barcodes[library_id][barcode] = sample_id

    for record in records:
        row_number = record["_row_number"]
        sample_id = str(record["sample_id"])
        is_control = bool(record["is_control"])
        control_id = record["control_id"]
        expected_motif = record["expected_motif"]
        if is_control:
            if record["assay_target"] != "IgG":
                raise ManifestValidationError(
                    f"row {row_number}: control {sample_id!r} must have assay_target IgG"
                )
            if control_id is not None:
                raise ManifestValidationError(
                    f"row {row_number}: control {sample_id!r} must not reference another control"
                )
            if expected_motif is not None:
                raise ManifestValidationError(
                    f"row {row_number}: expected_motif must be blank for controls"
                )
            continue

        if control_id is None:
            raise ManifestValidationError(
                f"row {row_number}: target {sample_id!r} is missing control_id"
            )
        if expected_motif is None:
            raise ManifestValidationError(
                f"row {row_number}: target {sample_id!r} is missing expected_motif"
            )
        control = samples.get(str(control_id))
        if control is None:
            raise ManifestValidationError(
                f"row {row_number}: target {sample_id!r} references missing control {control_id!r}"
            )
        if not bool(control["is_control"]):
            raise ManifestValidationError(
                f"row {row_number}: control {control_id!r} is not marked is_control"
            )
        if control["assay_target"] != "IgG":
            raise ManifestValidationError(
                f"row {row_number}: control {control_id!r} must have assay_target IgG"
            )
        if record["input_group"] != control["input_group"]:
            raise ManifestValidationError(
                f"row {row_number}: target and control must share input_group"
            )


def load_and_validate(path: Path, *, check_files: bool = True) -> list[dict[str, object]]:
    """Load a CSV manifest, validate it, and return normalized records."""

    manifest_path = Path(path).resolve()
    records = _read_rows(manifest_path, check_files)
    _validate_relationships(records)
    for record in records:
        del record["_row_number"]
    return records


def _write_json(records: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent
    )
    try:
        with os.fdopen(file_descriptor, "w", encoding="utf-8") as handle:
            json.dump(records, handle, indent=2)
            handle.write("\n")
        Path(temporary_name).replace(output)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate_parser = subparsers.add_parser("validate", help="validate and normalize a manifest")
    validate_parser.add_argument("--input", type=Path, required=True)
    validate_parser.add_argument("--output", type=Path, required=True)
    validate_parser.add_argument(
        "--skip-file-checks",
        action="store_true",
        help="validate manifest schema and relationships without requiring FASTQ files",
    )
    arguments = parser.parse_args(argv)

    try:
        records = load_and_validate(
            arguments.input, check_files=not arguments.skip_file_checks
        )
        _write_json(records, arguments.output)
    except ManifestValidationError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
