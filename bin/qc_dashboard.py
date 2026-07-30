#!/usr/bin/env python3
"""Parse and join the QC inputs used by the consolidated dashboard."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import html
import json
import math
import os
import re
import shutil
import statistics
import tempfile
from pathlib import Path
from typing import Mapping, Sequence

from motif_qc import read_ame


SCHEMA_VERSION = 1
GENERATOR_VERSION = "1.0.0"
TSS_BEFORE_BP = 3000
TSS_BIN_SIZE = 10
TSS_FLANK_BP = 100
TOP_MOTIF_LIMIT = 10
FAMILY_NAMES = (
    "demultiplex", "library", "insert_size", "peak",
    "peak_width", "tss", "motif", "ame",
)
FAMILY_STATUSES = (
    "computed", "skipped", "empty", "missing", "failed", "not_applicable",
)

QC_SUMMARY_COLUMNS = [
    "sample_id", "library_id", "input_group", "assay_target", "is_control",
    "control_id", "expected_motif", "total_read_pairs", "assigned_read_pairs",
    "ambiguous_read_pairs", "unassigned_read_pairs", "assigned_fraction",
    "ambiguous_fraction", "unassigned_fraction", "sample_assigned_reads",
    "sample_assignment_fraction", "raw_total_reads", "mapped_percent",
    "properly_paired_percent", "mapq_filtered_reads", "mapq_filtered_fragments",
    "mapq_filtered_fraction", "markdup_examined_reads", "duplicate_total",
    "duplicate_percent", "mitochondrial_percent", "estimated_library_size",
    "insert_size_total_pairs", "insert_size_min", "insert_size_q25",
    "insert_size_mean", "insert_size_median", "insert_size_q75",
    "insert_size_max", "peak_count", "total_covered_bases", "total_fragments",
    "fragments_in_peaks", "frip", "peak_width_min", "peak_width_q25",
    "peak_width_mean", "peak_width_median", "peak_width_q75", "peak_width_max",
    "tss_status", "tss_enrichment",
    "expected_motif_status", "best_motif_id", "best_adjusted_p_value",
    "ame_status", "warning_count",
]

TOP_MOTIF_COLUMNS = [
    "sample_id", "assay_target", "expected_motif", "rank", "motif_id",
    "motif_alt_id", "adjusted_p_value", "p_value", "effect",
    "positive_sequences",
]
TSS_PROFILE_COLUMNS = ["sample_id", "position_bp", "signal"]
DASHBOARD_OUTPUT_FILENAMES = (
    "qc_dashboard.html",
    "qc_summary.tsv",
    "qc_summary.json",
    "top_motifs.tsv",
    "tss_profiles.tsv",
)

DEMULTIPLEX_FIELDS = tuple(QC_SUMMARY_COLUMNS[7:16])
LIBRARY_FIELDS = (
    "raw_total_reads", "mapped_percent", "properly_paired_percent",
    "mapq_filtered_reads", "mapq_filtered_fragments", "mapq_filtered_fraction",
    "markdup_examined_reads", "duplicate_total", "duplicate_percent",
    "mitochondrial_percent", "estimated_library_size", "insert_size_total_pairs",
    "insert_size_min", "insert_size_q25", "insert_size_mean",
    "insert_size_median", "insert_size_q75", "insert_size_max",
)
LIBRARY_REQUIRED_VALUE_FIELDS = (
    "raw_total_reads", "mapped_percent", "properly_paired_percent",
    "mapq_filtered_reads", "mapq_filtered_fragments", "mapq_filtered_fraction",
    "markdup_examined_reads", "duplicate_total", "duplicate_percent",
    "mitochondrial_percent", "insert_size_total_pairs",
)
PEAK_PRODUCER_FIELDS = (
    "peak_count", "total_covered_bases",
    "peak_width_min", "peak_width_mean", "peak_width_median",
    "peak_width_max", "peak_width_q25", "peak_width_q75",
    "peak_score_count", "peak_score_min", "peak_score_q25",
    "peak_score_mean", "peak_score_median", "peak_score_q75",
    "peak_score_max", "signal_value_count", "signal_value_min",
    "signal_value_q25", "signal_value_mean", "signal_value_median",
    "signal_value_q75", "signal_value_max", "total_fragments",
    "fragments_in_peaks", "frip",
)
PEAK_REQUIRED_VALUE_FIELDS = (
    "peak_count", "total_covered_bases", "peak_score_count",
    "signal_value_count", "total_fragments", "fragments_in_peaks", "frip",
)
PEAK_FIELDS = (
    "peak_count", "total_covered_bases", "width_count", "width_min",
    "width_q25", "width_mean", "width_median", "width_q75", "width_max",
    "peak_score_count", "peak_score_min", "peak_score_q25", "peak_score_mean",
    "peak_score_median", "peak_score_q75", "peak_score_max",
    "signal_value_count", "signal_value_min", "signal_value_q25",
    "signal_value_mean", "signal_value_median", "signal_value_q75",
    "signal_value_max", "total_fragments", "fragments_in_peaks", "frip",
)


class DashboardInputError(ValueError):
    """Raised when a QC input cannot be joined without ambiguity."""


def finite_number(value: object, *, label: str, minimum: float | None = None) -> float:
    if isinstance(value, bool):
        raise DashboardInputError(f"{label} must be numeric")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise DashboardInputError(f"{label} must be numeric") from error
    if not math.isfinite(parsed) or (minimum is not None and parsed < minimum):
        raise DashboardInputError(f"{label} must be finite and >= {minimum}")
    return parsed


def _required_text(record: Mapping[str, object], field: str, *, label: str) -> str:
    value = record.get(field)
    if not isinstance(value, str) or not value.strip():
        raise DashboardInputError(f"{label} {field} is required")
    return value.strip()


def _is_blank(value: object) -> bool:
    return value is None or (isinstance(value, str) and not value.strip())


def _read_json(path: Path, *, label: str) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise DashboardInputError(f"cannot read {label} file {path}: {error}") from error


def load_metadata(path: Path) -> dict[str, dict[str, object]]:
    """Load uniquely keyed derived-sample metadata from a JSON array."""
    payload = _read_json(path, label="metadata")
    if not isinstance(payload, list):
        raise DashboardInputError("metadata must contain a JSON array")
    records: dict[str, dict[str, object]] = {}
    for index, raw_record in enumerate(payload, start=1):
        if not isinstance(raw_record, dict):
            raise DashboardInputError(f"metadata record {index} must be an object")
        sample_id = _required_text(raw_record, "sample_id", label="metadata")
        if sample_id in records:
            raise DashboardInputError(f"duplicate sample_id {sample_id}")
        record = {
            "sample_id": sample_id,
            "library_id": _required_text(raw_record, "library_id", label="metadata"),
            "assay_target": _required_text(raw_record, "assay_target", label="metadata"),
            "is_control": raw_record.get("is_control"),
            "input_group": raw_record.get("input_group"),
            "control_id": raw_record.get("control_id"),
            "expected_motif": raw_record.get("expected_motif"),
        }
        if not isinstance(record["is_control"], bool):
            raise DashboardInputError("metadata is_control must be boolean")
        records[sample_id] = record
    for sample_id, record in records.items():
        if record["is_control"]:
            if record["assay_target"] != "IgG":
                raise DashboardInputError(f"control {sample_id} assay_target must be IgG")
            if not _is_blank(record["control_id"]):
                raise DashboardInputError(f"control {sample_id} control_id must be blank")
            if not _is_blank(record["expected_motif"]):
                raise DashboardInputError(f"control {sample_id} expected_motif must be blank")
            continue
        control_id = record["control_id"]
        expected_motif = record["expected_motif"]
        if not isinstance(control_id, str) or not control_id.strip():
            raise DashboardInputError(f"target {sample_id} control_id is required")
        if not isinstance(expected_motif, str) or not expected_motif.strip():
            raise DashboardInputError(f"target {sample_id} expected_motif is required")
        control_id = control_id.strip()
        control = records.get(control_id)
        if control is None or not control["is_control"]:
            raise DashboardInputError(f"target {sample_id} control_id {control_id} must reference a control")
        if record["input_group"] != control["input_group"]:
            raise DashboardInputError(
                f"target {sample_id} control_id {control_id} must share input_group"
            )
    return dict(sorted(records.items()))


def _count(value: object, *, label: str) -> int | float:
    parsed = finite_number(value, label=label, minimum=0)
    return int(parsed) if parsed.is_integer() else parsed


def read_demultiplex_metrics(
    paths: Sequence[Path], metadata: Mapping[str, Mapping[str, object]] | None = None
) -> dict[str, dict[str, object]]:
    """Read physical-library count summaries and derive their fractions.

    Legacy inputs carry a ``library_id`` and can be parsed alone.  The current
    producer format omits it, so callers must supply validated sample metadata
    to recover the unique physical-library identity from ``assignment_counts``.
    """
    parsed: dict[str, dict[str, object]] = {}
    for path in paths:
        raw = _read_json(path, label="demultiplex metrics")
        if not isinstance(raw, dict):
            raise DashboardInputError(f"{path}: demultiplex metrics must be an object")
        counts = {
            name: _count(raw.get(name), label=f"demultiplex {name}")
            for name in ("total_reads", "assigned_reads", "ambiguous_reads", "unassigned_reads")
        }
        raw_assignments = raw.get("assignment_counts")
        if not isinstance(raw_assignments, dict):
            raise DashboardInputError(f"{path}: assignment_counts must be an object")
        assignment_counts: dict[str, int | float] = {}
        for sample_id, count in raw_assignments.items():
            if not isinstance(sample_id, str) or not sample_id.strip():
                raise DashboardInputError(f"{path}: assignment_counts sample_id is required")
            if metadata is not None and sample_id not in metadata:
                raise DashboardInputError(
                    f"assignment_counts reference unknown sample_id {sample_id}"
                )
            assignment_counts[sample_id] = _count(
                count, label=f"assignment count for {sample_id}"
            )
        if sum(assignment_counts.values()) != counts["assigned_reads"]:
            raise DashboardInputError("assignment_counts must sum to assigned_reads")
        if sum(counts[name] for name in ("assigned_reads", "ambiguous_reads", "unassigned_reads")) != counts["total_reads"]:
            raise DashboardInputError(
                "assigned_reads, ambiguous_reads, and unassigned_reads must sum to total_reads"
            )
        inferred_library_ids = (
            {
                _required_text(metadata[sample_id], "library_id", label=f"metadata {sample_id}")
                for sample_id in assignment_counts
            }
            if metadata is not None else set()
        )
        explicit_library_id = raw.get("library_id")
        if explicit_library_id is not None:
            if not isinstance(explicit_library_id, str) or not explicit_library_id.strip():
                raise DashboardInputError("demultiplex metrics library_id is required when supplied")
            library_id = explicit_library_id.strip()
            if metadata is not None:
                if library_id not in {
                    _required_text(record, "library_id", label=f"metadata {sample_id}")
                    for sample_id, record in metadata.items()
                }:
                    raise DashboardInputError(
                        f"demultiplex metrics reference unknown library_id {library_id}"
                    )
                if inferred_library_ids and inferred_library_ids != {library_id}:
                    raise DashboardInputError(
                        f"assignment_counts do not match library_id {library_id}"
                    )
        else:
            if metadata is None:
                raise DashboardInputError(
                    "validated metadata is required to infer library_id from assignment_counts"
                )
            if len(inferred_library_ids) != 1:
                raise DashboardInputError(
                    "assignment_counts map to multiple library_id values"
                    if inferred_library_ids else "assignment_counts cannot determine library_id"
                )
            library_id = next(iter(inferred_library_ids))
        if metadata is not None:
            expected_sample_ids = {
                sample_id
                for sample_id, record in metadata.items()
                if _required_text(
                    record, "library_id", label=f"metadata {sample_id}"
                ) == library_id
            }
            if set(assignment_counts) != expected_sample_ids:
                missing = sorted(expected_sample_ids - set(assignment_counts))
                extra = sorted(set(assignment_counts) - expected_sample_ids)
                detail = (
                    f"; missing {missing[0]}" if missing
                    else f"; unexpected {extra[0]}"
                )
                raise DashboardInputError(
                    f"assignment_counts for {library_id} must exactly match "
                    f"metadata samples{detail}"
                )
        if library_id in parsed:
            raise DashboardInputError(f"duplicate library_id {library_id}")
        total = counts["total_reads"]
        denominator = total or 1
        parsed[library_id] = {
            "library_id": library_id,
            **counts,
            "assignment_counts": dict(sorted(assignment_counts.items())),
            "assigned_fraction": counts["assigned_reads"] / denominator,
            "ambiguous_fraction": counts["ambiguous_reads"] / denominator,
            "unassigned_fraction": counts["unassigned_reads"] / denominator,
        }
    return dict(sorted(parsed.items()))


def _read_tsv(path: Path, *, label: str) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            headers = reader.fieldnames
            if not headers or any(not header for header in headers):
                raise DashboardInputError(f"{path}: {label} must have named TSV columns")
            if len(set(headers)) != len(headers):
                raise DashboardInputError(f"{path}: duplicate {label} column names")
            return headers, list(reader)
    except OSError as error:
        raise DashboardInputError(f"cannot read {label} file {path}: {error}") from error


def _read_single_sample_tsv(
    paths: Sequence[Path],
    *,
    label: str,
    expected_fields: Sequence[str] | None = None,
    required_value_fields: Sequence[str] = (),
) -> dict[str, dict[str, object]]:
    parsed: dict[str, dict[str, object]] = {}
    for path in paths:
        headers, rows = _read_tsv(path, label=label)
        if (
            expected_fields is not None
            and headers != ["sample_id", *expected_fields]
        ):
            raise DashboardInputError(
                f"{path}: {label} columns must exactly match the producer schema"
            )
        if "sample_id" not in headers or len(rows) != 1:
            raise DashboardInputError(f"{path}: {label} must contain exactly one sample_id row")
        raw = rows[0]
        sample_id = _required_text(raw, "sample_id", label=label)
        if sample_id in parsed:
            raise DashboardInputError(f"duplicate sample_id {sample_id}")
        record: dict[str, object] = {"sample_id": sample_id}
        for name, value in raw.items():
            if name == "sample_id":
                continue
            if name in required_value_fields and _is_blank(value):
                raise DashboardInputError(
                    f"{path}: required {label} field {name} is empty"
                )
            record[name] = None if _is_blank(value) else finite_number(
                value, label=f"{sample_id} {name}", minimum=0
            )
        parsed[sample_id] = record
    return dict(sorted(parsed.items()))


def read_library_metrics(paths: Sequence[Path]) -> dict[str, dict[str, object]]:
    """Read one library-QC summary row for every derived sample."""
    return _read_single_sample_tsv(
        paths,
        label="library metrics",
        expected_fields=LIBRARY_FIELDS,
        required_value_fields=LIBRARY_REQUIRED_VALUE_FIELDS,
    )


def read_peak_metrics(paths: Sequence[Path]) -> dict[str, dict[str, object]]:
    """Convert PEAK_QC metric/value tables into sample-keyed records."""
    parsed: dict[str, dict[str, object]] = {}
    for path in paths:
        headers, rows = _read_tsv(path, label="peak metrics")
        if headers != ["metric", "value"]:
            raise DashboardInputError(f"{path}: peak metrics must have metric and value columns")
        raw_metrics: dict[str, str] = {}
        for row in rows:
            metric = row["metric"].strip()
            if not metric:
                raise DashboardInputError(f"{path}: peak metric is required")
            if metric in raw_metrics:
                raise DashboardInputError(f"duplicate peak metric {metric}")
            raw_metrics[metric] = row["value"]
        expected_metrics = {"sample_id", *PEAK_PRODUCER_FIELDS}
        missing_metrics = sorted(expected_metrics - set(raw_metrics))
        if missing_metrics:
            raise DashboardInputError(
                f"{path}: peak metrics are missing required producer fields: "
                f"{missing_metrics[0]}"
            )
        unexpected_metrics = sorted(set(raw_metrics) - expected_metrics)
        if unexpected_metrics:
            raise DashboardInputError(
                f"{path}: peak metrics contain unexpected producer field "
                f"{unexpected_metrics[0]}"
            )
        sample_id = raw_metrics.get("sample_id", "").strip()
        if not sample_id:
            raise DashboardInputError(f"{path}: peak sample_id is required")
        if sample_id in parsed:
            raise DashboardInputError(f"duplicate sample_id {sample_id}")
        record: dict[str, object] = {"sample_id": sample_id}
        for metric, value in raw_metrics.items():
            if metric != "sample_id":
                if metric in PEAK_REQUIRED_VALUE_FIELDS and not value.strip():
                    raise DashboardInputError(
                        f"{path}: required peak metric {metric} is empty"
                    )
                internal_metric = (
                    metric.removeprefix("peak_")
                    if metric.startswith("peak_width_")
                    else metric
                )
                record[internal_metric] = None if not value.strip() else finite_number(
                    value, label=f"{sample_id} peak {metric}", minimum=0
                )
        record["width_count"] = record["peak_count"]
        if record["peak_count"]:
            for metric in (
                "width_min", "width_mean", "width_median",
                "width_max", "width_q25", "width_q75",
            ):
                if record[metric] is None:
                    raise DashboardInputError(
                        f"{path}: peak metric {metric} is required when peaks exist"
                    )
        parsed[sample_id] = record
    return dict(sorted(parsed.items()))


def _sample_id_from_filename(path: Path, suffix: str, *, label: str) -> str:
    if not path.name.endswith(suffix):
        raise DashboardInputError(f"{path}: unexpected {label} filename")
    sample_id = path.name[:-len(suffix)]
    if not sample_id:
        raise DashboardInputError(f"{path}: {label} filename sample_id is required")
    return sample_id


def _nonnegative_integer(value: object, *, label: str) -> int:
    parsed = finite_number(value, label=label, minimum=0)
    if not parsed.is_integer():
        raise DashboardInputError(f"{label} must be an integer")
    return int(parsed)


def read_insert_size_distributions(
    paths: Sequence[Path],
) -> dict[str, list[dict[str, int]]]:
    """Read sample-keyed insert-size histograms, retaining header-only empties."""
    records: dict[str, list[dict[str, int]]] = {}
    suffix = ".insert_size_distribution.tsv"
    for path in paths:
        sample_id = _sample_id_from_filename(path, suffix, label="insert-size distribution")
        if sample_id in records:
            raise DashboardInputError(f"duplicate insert-size distribution sample_id {sample_id}")
        headers, rows = _read_tsv(path, label="insert-size distribution")
        if headers != ["sample_id", "insert_size", "pair_count"]:
            raise DashboardInputError(
                f"{path}: insert-size distribution must have sample_id, "
                "insert_size, and pair_count columns"
            )
        parsed_rows: list[dict[str, int]] = []
        seen_sizes: set[int] = set()
        for row_number, row in enumerate(rows, start=2):
            row_sample_id = _required_text(
                row, "sample_id", label=f"{path}:{row_number}"
            )
            if row_sample_id != sample_id:
                raise DashboardInputError(
                    f"{path}:{row_number}: distribution sample_id {row_sample_id} "
                    f"does not match filename sample_id {sample_id}"
                )
            insert_size = _nonnegative_integer(
                row["insert_size"], label=f"{path}:{row_number} insert_size"
            )
            pair_count = _nonnegative_integer(
                row["pair_count"], label=f"{path}:{row_number} pair_count"
            )
            if insert_size in seen_sizes:
                raise DashboardInputError(
                    f"{path}: duplicate insert_size {insert_size}"
                )
            seen_sizes.add(insert_size)
            parsed_rows.append({"insert_size": insert_size, "pair_count": pair_count})
        records[sample_id] = sorted(
            parsed_rows, key=lambda row: row["insert_size"]
        )
    return dict(sorted(records.items()))


def read_peak_width_distributions(
    paths: Sequence[Path],
) -> dict[str, list[dict[str, int]]]:
    """Read target peak-width histograms, retaining header-only empty results."""
    records: dict[str, list[dict[str, int]]] = {}
    suffix = ".peak_qc.width_histogram.tsv"
    for path in paths:
        sample_id = _sample_id_from_filename(path, suffix, label="peak-width distribution")
        if sample_id in records:
            raise DashboardInputError(f"duplicate peak-width distribution sample_id {sample_id}")
        headers, rows = _read_tsv(path, label="peak-width distribution")
        if headers != ["width", "peak_count"]:
            raise DashboardInputError(
                f"{path}: peak-width distribution must have width and peak_count columns"
            )
        parsed_rows: list[dict[str, int]] = []
        seen_widths: set[int] = set()
        for row_number, row in enumerate(rows, start=2):
            width = _nonnegative_integer(
                row["width"], label=f"{path}:{row_number} width"
            )
            peak_count = _nonnegative_integer(
                row["peak_count"], label=f"{path}:{row_number} peak_count"
            )
            if width in seen_widths:
                raise DashboardInputError(f"{path}: duplicate peak width {width}")
            seen_widths.add(width)
            parsed_rows.append({"width": width, "peak_count": peak_count})
        records[sample_id] = sorted(parsed_rows, key=lambda row: row["width"])
    return dict(sorted(records.items()))


def read_fragments_per_peak_distributions(
    paths: Sequence[Path],
) -> dict[str, list[dict[str, int]]]:
    """Compact strict per-peak producer rows into sample-keyed histograms."""
    records: dict[str, list[dict[str, int]]] = {}
    suffix = ".peak_qc.fragments_per_peak.tsv"
    expected_fields = [
        "chrom", "start", "end", "peak_name", "width", "score",
        "signal_value", "fragment_count",
    ]
    for path in paths:
        sample_id = _sample_id_from_filename(
            path, suffix, label="fragments-per-peak distribution"
        )
        if sample_id in records:
            raise DashboardInputError(
                f"duplicate fragments-per-peak distribution sample_id {sample_id}"
            )
        headers, rows = _read_tsv(path, label="fragments-per-peak distribution")
        if headers != expected_fields:
            raise DashboardInputError(
                f"{path}: fragments-per-peak distribution columns must exactly "
                "match the producer schema"
            )
        counts = Counter(
            _nonnegative_integer(
                row["fragment_count"],
                label=f"{path}:{row_number} fragment_count",
            )
            for row_number, row in enumerate(rows, start=2)
        )
        records[sample_id] = [
            {"fragment_count": fragment_count, "peak_count": counts[fragment_count]}
            for fragment_count in sorted(counts)
        ]
    return dict(sorted(records.items()))


def read_tss_profile(
    path: Path, *, before_bp: int = TSS_BEFORE_BP, bin_size: int = TSS_BIN_SIZE
) -> list[tuple[int, float | None]]:
    """Read a pinned deepTools 3.5.5 ``plotProfile --outFileNameData`` table."""
    if before_bp <= 0 or bin_size <= 0 or (2 * before_bp) % bin_size:
        raise DashboardInputError("TSS profile window must divide evenly into bins")
    expected_bins = (2 * before_bp) // bin_size
    try:
        rows = [
            line.split("\t")
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    except OSError as error:
        raise DashboardInputError(f"cannot read TSS profile {path}: {error}") from error
    if len(rows) != 3:
        raise DashboardInputError(
            f"{path}: TSS profile must contain exactly three non-comment rows"
        )
    expected_fields = expected_bins + 2
    for line_number, row in enumerate(rows, start=1):
        if len(row) != expected_fields:
            raise DashboardInputError(
                f"{path}:{line_number}: TSS profile row must contain "
                f"{expected_fields} tab-separated fields"
            )
    labels, bins, data_row = rows
    if labels[:2] != ["bin labels", ""]:
        raise DashboardInputError(
            f"{path}: first TSS profile row must start with 'bin labels' and a blank field"
        )
    if bins[:2] != ["bins", ""]:
        raise DashboardInputError(
            f"{path}: second TSS profile row must start with 'bins' and a blank field"
        )
    try:
        numeric_bins = [float(value) for value in bins[2:]]
    except ValueError as error:
        raise DashboardInputError(f"{path}: TSS profile bins must be integers") from error
    if any(
        not math.isfinite(value) or not value.is_integer() for value in numeric_bins
    ):
        raise DashboardInputError(f"{path}: TSS profile bins must be integers")
    parsed_bins = [int(value) for value in numeric_bins]
    if parsed_bins != list(range(1, expected_bins + 1)):
        raise DashboardInputError(
            f"{path}: TSS profile bins must be exactly 1 through {expected_bins}"
        )
    if not data_row[0].strip() or not data_row[1].strip():
        raise DashboardInputError(
            f"{path}: TSS profile sample and group labels are required"
        )
    values: list[float | None] = []
    for index, field in enumerate(data_row[2:], start=1):
        if field.strip() == "--":
            values.append(None)
            continue
        try:
            value = float(field)
        except ValueError as error:
            raise DashboardInputError(
                f"{path}: TSS profile bin {index} must be numeric or missing"
            ) from error
        values.append(value if math.isfinite(value) else None)
    return [
        (-before_bp + index * bin_size, value)
        for index, value in enumerate(values)
    ]


def _tss_enrichment_result(
    profile: Sequence[tuple[int, float | None]], *, flank_bp: int = TSS_FLANK_BP
) -> tuple[float | None, str]:
    """Return the scalar and a diagnostic suitable for structured warnings."""
    if flank_bp <= 0 or flank_bp % TSS_BIN_SIZE:
        raise DashboardInputError("TSS flank_bp must be a positive multiple of bin size")
    flank_bins = flank_bp // TSS_BIN_SIZE
    if len(profile) < 2 * flank_bins + 1:
        raise DashboardInputError("TSS profile is too short for the requested flanks")
    center_raw = profile[len(profile) // 2][1]
    if center_raw is None or not math.isfinite(float(center_raw)):
        return None, "non_finite_center"
    flank_raw = [value for _, value in profile[:flank_bins]]
    flank_raw += [value for _, value in profile[-flank_bins:]]
    finite_flanks = [
        float(value)
        for value in flank_raw
        if value is not None and math.isfinite(float(value))
    ]
    if not finite_flanks:
        return None, "non_finite_flank"
    flank_mean = statistics.fmean(finite_flanks)
    if flank_mean == 0:
        return None, "zero_flank"
    score = float(center_raw) / flank_mean
    if not math.isfinite(score):
        return None, "non_finite_ratio"
    if len(finite_flanks) != len(flank_raw):
        return score, "partial_non_finite_flank"
    return score, "computed"


def calculate_tss_enrichment(
    profile: Sequence[tuple[int, float | None]], *, flank_bp: int = TSS_FLANK_BP
) -> float | None:
    """Return center signal divided by the two terminal TSS flank means."""
    return _tss_enrichment_result(profile, flank_bp=flank_bp)[0]


def _ame_record_dict(record: object) -> dict[str, object]:
    """Serialize one parsed AME record for dashboard data structures."""
    return {
        "motif_id": record.motif_id,
        "motif_alt_id": record.motif_alt_id,
        "adjusted_p_value": record.adjusted_p_value,
        "p_value": record.p_value,
        "effect": record.effect,
        "positive_sequences": record.positive_sequences,
        "rank": record.rank,
    }


def read_all_ame(path: Path) -> list[dict[str, object]]:
    """Return all usable AME records in deterministic significance order."""
    try:
        records = read_ame(path)
    except ValueError as error:
        raise DashboardInputError(str(error)) from error
    usable = [record for record in records if record.motif_id != "__NO_PEAKS__"]
    usable.sort(key=lambda record: (
        record.adjusted_p_value,
        record.rank if record.rank is not None else math.inf,
        record.motif_id,
        record.motif_alt_id,
    ))
    return [_ame_record_dict(record) for record in usable]


def read_top_ame(path: Path, *, limit: int = TOP_MOTIF_LIMIT) -> list[dict[str, object]]:
    """Return up to ``limit`` AME records in deterministic significance order."""
    if limit < 0:
        raise DashboardInputError("AME motif limit must be non-negative")
    return read_all_ame(path)[:limit]


def motif_tokens(*values: object) -> frozenset[str]:
    """Return punctuation-delimited, case-insensitive motif tokens."""
    return frozenset(
        token
        for value in values
        if value is not None
        for token in re.findall(r"[A-Za-z0-9]+", str(value).casefold())
    )


def is_cognate_motif(
    expected_motif: object, motif_id: object, motif_alt_id: object
) -> bool:
    """Return whether the expected TF occurs as a complete motif token."""
    expected_tokens = motif_tokens(expected_motif)
    return bool(expected_tokens & motif_tokens(motif_id, motif_alt_id))


def build_motif_heatmap(
    samples: Sequence[Mapping[str, object]],
    *,
    noncognate_limit: int = 15,
    significance_cap: float = 60.0,
) -> dict[str, object]:
    """Build a deterministic target-only AME significance matrix."""
    targets = sorted(
        (sample for sample in samples if not bool(sample.get("is_control"))),
        key=lambda sample: str(sample.get("sample_id", "")),
    )
    sample_ids = [str(sample.get("sample_id", "")) for sample in targets]
    expected_motifs = sorted(
        {
            str(sample.get("expected_motif")).strip()
            for sample in targets
            if motif_tokens(sample.get("expected_motif"))
        },
        key=lambda value: (value.casefold(), value),
    )
    records_by_sample: dict[
        str, dict[tuple[str, str], Mapping[str, object]]
    ] = {}
    key_records: dict[tuple[str, str], list[Mapping[str, object]]] = {}
    display_candidates: dict[tuple[str, str], set[tuple[str, str]]] = {}
    cognate_groups: dict[tuple[str, str], str] = {}

    def display_key(record: Mapping[str, object]) -> tuple[str, str]:
        return (
            str(record.get("motif_id") or "").strip(),
            str(record.get("motif_alt_id") or "").strip(),
        )

    def normalized_key(key: tuple[str, str]) -> tuple[str, str]:
        return key[0].casefold(), key[1].casefold()

    def record_priority(record: Mapping[str, object]) -> tuple[object, ...]:
        adjusted = record.get("adjusted_p_value")
        adjusted_sort = (
            float(adjusted)
            if isinstance(adjusted, (int, float))
            and math.isfinite(float(adjusted))
            else math.inf
        )
        rank = record.get("rank")
        rank_sort = (
            float(rank)
            if isinstance(rank, (int, float)) and math.isfinite(float(rank))
            else math.inf
        )
        key = display_key(record)
        return adjusted_sort, rank_sort, key[1], key[0]

    for sample in targets:
        sample_id = str(sample.get("sample_id", ""))
        sample_records: dict[tuple[str, str], Mapping[str, object]] = {}
        raw_records = sample.get("ame_motifs", [])
        if not isinstance(raw_records, Sequence) or isinstance(
            raw_records, (str, bytes)
        ):
            raw_records = []
        for record in raw_records:
            if not isinstance(record, Mapping):
                continue
            display = display_key(record)
            key = normalized_key(display)
            current_record = sample_records.get(key)
            if (
                current_record is None
                or record_priority(record) < record_priority(current_record)
            ):
                sample_records[key] = record
            key_records.setdefault(key, []).append(record)
            display_candidates.setdefault(key, set()).add(display)
        records_by_sample[sample_id] = sample_records

    display_by_key = {
        key: min(
            candidates,
            key=lambda display: (
                display[1].casefold(),
                display[0].casefold(),
                display[1],
                display[0],
            ),
        )
        for key, candidates in display_candidates.items()
    }
    for key, display in display_by_key.items():
        matching_groups = [
            expected.casefold()
            for expected in expected_motifs
            if is_cognate_motif(expected, *display)
        ]
        if matching_groups:
            cognate_groups[key] = min(matching_groups)

    forced_normalized_keys = sorted(
        cognate_groups,
        key=lambda key: (
            cognate_groups[key],
            display_by_key[key][1].casefold(),
            display_by_key[key][0].casefold(),
            display_by_key[key],
        ),
    )

    def best_adjusted(key: tuple[str, str]) -> float:
        values = []
        for record in key_records[key]:
            value = record.get("adjusted_p_value")
            if isinstance(value, (int, float)) and math.isfinite(float(value)):
                values.append(float(value))
        return min(values) if values else math.inf

    noncognate_normalized_keys = sorted(
        (key for key in key_records if key not in cognate_groups),
        key=lambda key: (
            best_adjusted(key),
            display_by_key[key][1].casefold(),
            display_by_key[key][0].casefold(),
            display_by_key[key],
        ),
    )[:noncognate_limit]
    normalized_keys = [
        *forced_normalized_keys, *noncognate_normalized_keys,
    ]
    forced_cognate_keys = [
        display_by_key[key] for key in forced_normalized_keys
    ]
    noncognate_keys = [
        display_by_key[key] for key in noncognate_normalized_keys
    ]
    motif_keys = [display_by_key[key] for key in normalized_keys]
    cells: dict[tuple[str, tuple[str, str]], dict[str, object]] = {}
    for sample in targets:
        sample_id = str(sample.get("sample_id", ""))
        expected = sample.get("expected_motif")
        for key in normalized_keys:
            display = display_by_key[key]
            record = records_by_sample[sample_id].get(key)
            adjusted = (
                record.get("adjusted_p_value")
                if isinstance(record, Mapping) else None
            )
            significant = (
                isinstance(adjusted, (int, float))
                and math.isfinite(float(adjusted))
                and 0 <= float(adjusted) <= 0.1
            )
            if significant and float(adjusted) == 0:
                score = significance_cap
                label = f">{format(significance_cap, 'g')}"
            elif significant:
                score = min(-math.log10(float(adjusted)), significance_cap)
                label = format(score, ".3g")
            else:
                score = 0.0
                label = "ns"
            cells[(sample_id, display)] = {
                "score": score,
                "label": label,
                "adjusted_p_value": adjusted,
                "outlined": is_cognate_motif(expected, *display),
            }

    return {
        "sample_ids": sample_ids,
        "motif_keys": motif_keys,
        "forced_cognate_keys": forced_cognate_keys,
        "noncognate_keys": noncognate_keys,
        "cells": cells,
    }


def _family_state(status: str, reason: str | None = None) -> dict[str, object]:
    if status not in FAMILY_STATUSES:
        raise DashboardInputError(f"unsupported QC family status {status}")
    return {"status": status, "reason": reason}


def _empty_sample(record: Mapping[str, object]) -> dict[str, object]:
    return {
        **dict(record),
        "sample_kind": "control" if record["is_control"] else "target",
        "demultiplex": {field: None for field in DEMULTIPLEX_FIELDS},
        "library": {"insert_size_distribution": None},
        "peak": {
            "frip": None,
            "width_distribution": None,
            "fragments_per_peak_distribution": None,
        },
        "tss": {
            "status": None, "enrichment": None, "profile": None,
            "score_status": None,
        },
        "motif": {
            "status": "not_applicable_control" if record["is_control"] else None,
            "best_motif_id": None,
            "best_adjusted_p_value": None,
            "ame_status": "not_applicable_control" if record["is_control"] else None,
        },
        "availability": {
            family: _family_state("missing", "not_evaluated")
            for family in FAMILY_NAMES
        },
        "ame_motifs": [],
        "top_motifs": [],
        "warnings": [],
    }


def _set_availability(
    sample: dict[str, object], family: str, status: str, reason: str | None = None
) -> None:
    availability = sample["availability"]
    if not isinstance(availability, dict):
        raise DashboardInputError("sample availability must be an object")
    availability[family] = _family_state(status, reason)


def _add_warning(
    sample: dict[str, object], family: str, status: str, message: str
) -> None:
    warnings = sample["warnings"]
    if not isinstance(warnings, list):
        raise DashboardInputError("sample warnings must be a list")
    warnings.append({
        "sample_id": sample["sample_id"],
        "family": family,
        "status": status,
        "message": message,
    })


def _status_from_optional_result(status: object) -> tuple[str, str | None]:
    if status in {"computed", "pass", "not_significant", "motif_not_found"}:
        return "computed", None
    if status in {"empty", "no_peaks"}:
        return "empty", "no_usable_results"
    if status in {"skipped", "skipped_no_annotation", "skipped_no_database"}:
        return "skipped", str(status)
    if status == "failed":
        return "failed", "producer_failed"
    if status is None:
        return "missing", "status_missing"
    return "failed", f"unsupported_status_{status}"


def _aggregate_availability(
    samples: Sequence[Mapping[str, object]],
) -> dict[str, dict[str, object]]:
    aggregates: dict[str, dict[str, object]] = {}
    priority = ("failed", "missing", "empty", "skipped", "computed", "not_applicable")
    for family in FAMILY_NAMES:
        counts = {status: 0 for status in FAMILY_STATUSES}
        for sample in samples:
            availability = sample.get("availability", {})
            entry = availability.get(family, {}) if isinstance(availability, Mapping) else {}
            status = entry.get("status") if isinstance(entry, Mapping) else "missing"
            if status not in counts:
                status = "failed"
            counts[str(status)] += 1
        run_status = next(
            status for status in priority if counts[status]
        ) if samples else "missing"
        if counts["computed"] and not any(
            counts[status] for status in ("failed", "missing", "empty", "skipped")
        ):
            run_status = "computed"
        aggregates[family] = {"status": run_status, **counts}
    return aggregates


def build_report_data(
    metadata: Mapping[str, Mapping[str, object]],
    demultiplex: Mapping[str, Mapping[str, object]],
    libraries: Mapping[str, Mapping[str, object]],
    peaks: Mapping[str, Mapping[str, object]],
    tss: Mapping[str, Mapping[str, object]],
    motifs: Mapping[str, Mapping[str, object]],
    top_motifs: Mapping[str, Sequence[Mapping[str, object]]],
    annotation_status: str,
    *,
    insert_sizes: Mapping[str, Sequence[Mapping[str, object]]] | None = None,
    peak_widths: Mapping[str, Sequence[Mapping[str, object]]] | None = None,
    fragments_per_peak: (
        Mapping[str, Sequence[Mapping[str, object]]] | None
    ) = None,
    ame_motifs: Mapping[str, Sequence[Mapping[str, object]]] | None = None,
    motif_analysis_status: str = "auto",
) -> dict[str, object]:
    """Initialize optional QC families, then overlay validated sample metrics."""
    if annotation_status not in {
        "skipped_no_annotation", "computed_bed", "computed_gtf",
    }:
        raise DashboardInputError(f"unsupported annotation status {annotation_status}")
    if motif_analysis_status == "auto":
        motif_analysis_status = (
            "computed" if motifs or top_motifs or ame_motifs else "skipped_no_database"
        )
    if motif_analysis_status not in {"computed", "skipped_no_database"}:
        raise DashboardInputError(
            f"unsupported motif analysis status {motif_analysis_status}"
        )
    insert_sizes = {} if insert_sizes is None else insert_sizes
    peak_widths = {} if peak_widths is None else peak_widths
    fragments_per_peak = {} if fragments_per_peak is None else fragments_per_peak
    ame_motifs = {} if ame_motifs is None else ame_motifs
    samples_by_id = {sample_id: _empty_sample(record) for sample_id, record in metadata.items()}
    sample_ids = set(samples_by_id)
    for family_name, family in (("library", libraries), ("peak", peaks), ("tss", tss),
                                ("motif", motifs), ("top motifs", top_motifs),
                                ("insert-size", insert_sizes),
                                ("peak-width", peak_widths),
                                ("fragments-per-peak", fragments_per_peak),
                                ("AME motifs", ame_motifs)):
        unknown = sorted(set(family) - sample_ids)
        if unknown:
            raise DashboardInputError(f"{family_name} metrics reference unknown sample_id {unknown[0]}")
    library_ids = {sample["library_id"] for sample in samples_by_id.values()}
    for library_id, record in demultiplex.items():
        if library_id not in library_ids:
            raise DashboardInputError(f"demultiplex metrics reference unknown library_id {library_id}")
        counts = record.get("assignment_counts", {})
        if not isinstance(counts, Mapping):
            raise DashboardInputError(f"demultiplex {library_id} assignment_counts must be an object")
        unknown = sorted(set(counts) - sample_ids)
        if unknown:
            raise DashboardInputError(
                f"demultiplex metrics reference unknown sample_id {unknown[0]}"
            )
    for sample_id, sample in samples_by_id.items():
        library_id = sample["library_id"]
        demux = demultiplex.get(library_id)
        if demux is None:
            _set_availability(sample, "demultiplex", "missing", "library_metrics_missing")
            _add_warning(
                sample, "demultiplex", "missing",
                f"{sample_id}: missing demultiplex metrics for {library_id}",
            )
        else:
            counts = demux["assignment_counts"]
            total_assigned = demux.get("assigned_reads")
            sample_count = counts.get(sample_id)
            sample["demultiplex"].update({
                "total_read_pairs": demux.get("total_reads"),
                "assigned_read_pairs": total_assigned,
                "ambiguous_read_pairs": demux.get("ambiguous_reads"),
                "unassigned_read_pairs": demux.get("unassigned_reads"),
                "assigned_fraction": demux.get("assigned_fraction"),
                "ambiguous_fraction": demux.get("ambiguous_fraction"),
                "unassigned_fraction": demux.get("unassigned_fraction"),
                "sample_assigned_reads": sample_count,
                "sample_assignment_fraction": (
                    sample_count / total_assigned if sample_count is not None and total_assigned else None
                ),
            })
            _set_availability(sample, "demultiplex", "computed")
        if sample_id in libraries:
            sample["library"].update(libraries[sample_id])
            _set_availability(sample, "library", "computed")
        else:
            _set_availability(sample, "library", "missing", "metrics_missing")
            _add_warning(
                sample, "library", "missing",
                f"{sample_id}: missing library metrics",
            )

        if sample_id in insert_sizes:
            distribution = [dict(row) for row in insert_sizes[sample_id]]
            sample["library"]["insert_size_distribution"] = distribution
            insert_status = "computed" if distribution else "empty"
            _set_availability(
                sample, "insert_size", insert_status,
                None if distribution else "distribution_empty",
            )
            if not distribution:
                _add_warning(
                    sample, "insert_size", "empty",
                    f"{sample_id}: insert-size distribution is empty",
                )
        else:
            _set_availability(sample, "insert_size", "missing", "distribution_missing")
            _add_warning(
                sample, "insert_size", "missing",
                f"{sample_id}: missing insert-size distribution",
            )

        if sample["is_control"]:
            _set_availability(sample, "peak", "not_applicable", "control_library")
            _set_availability(sample, "peak_width", "not_applicable", "control_library")
        else:
            if sample_id in peaks:
                sample["peak"].update(peaks[sample_id])
                peak_status = "empty" if peaks[sample_id].get("peak_count") == 0 else "computed"
                _set_availability(
                    sample, "peak", peak_status,
                    "no_peaks" if peak_status == "empty" else None,
                )
                if peak_status == "empty":
                    _add_warning(
                        sample, "peak", "empty",
                        f"{sample_id}: peak result is empty",
                    )
            else:
                _set_availability(sample, "peak", "missing", "metrics_missing")
                _add_warning(
                    sample, "peak", "missing",
                    f"{sample_id}: missing peak metrics",
                )
            if sample_id in peak_widths:
                distribution = [dict(row) for row in peak_widths[sample_id]]
                sample["peak"]["width_distribution"] = distribution
                width_status = "computed" if distribution else "empty"
                _set_availability(
                    sample, "peak_width", width_status,
                    None if distribution else "distribution_empty",
                )
                if not distribution:
                    _add_warning(
                        sample, "peak_width", "empty",
                        f"{sample_id}: peak-width distribution is empty",
                    )
            else:
                _set_availability(
                    sample, "peak_width", "missing", "distribution_missing"
                )
                _add_warning(
                    sample, "peak_width", "missing",
                    f"{sample_id}: missing peak-width distribution",
                )
            if sample_id in fragments_per_peak:
                sample["peak"]["fragments_per_peak_distribution"] = [
                    dict(row) for row in fragments_per_peak[sample_id]
                ]

        if annotation_status == "skipped_no_annotation":
            sample["tss"]["status"] = "skipped_no_annotation"
            _set_availability(sample, "tss", "skipped", "annotation_not_provided")
            _add_warning(
                sample, "tss", "skipped",
                f"{sample_id}: TSS analysis skipped because no annotation was provided",
            )
        elif sample_id not in tss:
            _set_availability(sample, "tss", "missing", "metrics_missing")
            _add_warning(
                sample, "tss", "missing",
                f"{sample_id}: missing TSS metrics",
            )
        else:
            sample["tss"].update(tss[sample_id])
            tss_status, tss_reason = _status_from_optional_result(
                sample["tss"].get("status")
            )
            if tss_status == "computed" and sample["tss"].get("profile") is None:
                tss_status, tss_reason = "failed", "profile_missing"
                sample["tss"]["status"] = "failed"
                sample["tss"]["reason"] = tss_reason
            _set_availability(sample, "tss", tss_status, tss_reason)
            if tss_status != "computed":
                _add_warning(
                    sample, "tss", tss_status,
                    f"{sample_id}: TSS result is {tss_status}"
                    + (f" ({tss_reason})" if tss_reason else ""),
                )
            score_status = sample["tss"].get("score_status")
            if tss_status == "computed" and score_status not in {None, "computed"}:
                _add_warning(
                    sample, "tss", "computed",
                    f"{sample_id}: TSS enrichment is NA or partial ({score_status})",
                )

        if sample["is_control"]:
            _set_availability(sample, "motif", "not_applicable", "control_library")
            _set_availability(sample, "ame", "not_applicable", "control_library")
        elif motif_analysis_status == "skipped_no_database":
            sample["motif"]["status"] = "skipped_no_database"
            sample["motif"]["ame_status"] = "skipped_no_database"
            _set_availability(sample, "motif", "skipped", "motif_database_not_provided")
            _set_availability(sample, "ame", "skipped", "motif_database_not_provided")
            _add_warning(
                sample, "motif", "skipped",
                f"{sample_id}: motif analysis skipped because no motif database was provided",
            )
            _add_warning(
                sample, "ame", "skipped",
                f"{sample_id}: AME analysis skipped because no motif database was provided",
            )
        else:
            if sample_id in motifs:
                sample["motif"].update(motifs[sample_id])
            motif_status, motif_reason = _status_from_optional_result(
                sample["motif"].get("status")
            )
            _set_availability(sample, "motif", motif_status, motif_reason)
            if motif_status != "computed":
                _add_warning(
                    sample, "motif", motif_status,
                    f"{sample_id}: expected-motif result is {motif_status}"
                    + (f" ({motif_reason})" if motif_reason else ""),
                )
            ame_status, ame_reason = _status_from_optional_result(
                sample["motif"].get("ame_status")
            )
            _set_availability(sample, "ame", ame_status, ame_reason)
            if ame_status != "computed":
                _add_warning(
                    sample, "ame", ame_status,
                    f"{sample_id}: AME result is {ame_status}"
                    + (f" ({ame_reason})" if ame_reason else ""),
                )
            sample["top_motifs"] = [
                dict(item) for item in top_motifs.get(sample_id, ())
            ]
            sample["ame_motifs"] = [
                dict(item) for item in ame_motifs.get(sample_id, ())
            ]
    samples = [samples_by_id[sample_id] for sample_id in sorted(samples_by_id)]
    warnings = [warning for sample in samples for warning in sample["warnings"]]
    target_count = sum(not bool(sample["is_control"]) for sample in samples)
    control_count = len(samples) - target_count
    return {
        "schema_version": SCHEMA_VERSION,
        "generator_version": GENERATOR_VERSION,
        "annotation_status": annotation_status,
        "motif_analysis_status": motif_analysis_status,
        "counts": {
            "samples": len(samples),
            "targets": target_count,
            "controls": control_count,
            "warnings": len(warnings),
        },
        "availability": _aggregate_availability(samples),
        "samples": samples,
        "samples_by_id": samples_by_id,
        "warnings": warnings,
    }


def format_value(value: object) -> str:
    """Format a report value for people; machine TSV writers retain empty nulls."""
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return format(value, ".6g")
    return str(value)


def _samples(data: Mapping[str, object]) -> list[Mapping[str, object]]:
    samples = data.get("samples", [])
    if not isinstance(samples, list):
        raise DashboardInputError("report data samples must be a list")
    return sorted(
        (sample for sample in samples if isinstance(sample, Mapping)),
        key=lambda sample: str(sample.get("sample_id", "")),
    )


def _nested(sample: Mapping[str, object], key: str) -> Mapping[str, object]:
    value = sample.get(key, {})
    return value if isinstance(value, Mapping) else {}


def _summary_row(sample: Mapping[str, object]) -> dict[str, object]:
    demultiplex = _nested(sample, "demultiplex")
    library = _nested(sample, "library")
    peak = _nested(sample, "peak")
    tss = _nested(sample, "tss")
    motif = _nested(sample, "motif")
    row = {key: sample.get(key) for key in (
        "sample_id", "library_id", "input_group", "assay_target", "is_control",
        "control_id", "expected_motif",
    )}
    row.update({key: demultiplex.get(key) for key in QC_SUMMARY_COLUMNS if key in demultiplex})
    row.update({key: library.get(key) for key in QC_SUMMARY_COLUMNS if key in library})
    row.update({key: peak.get(key) for key in QC_SUMMARY_COLUMNS if key in peak})
    for suffix in ("min", "q25", "mean", "median", "q75", "max"):
        row[f"peak_width_{suffix}"] = peak.get(f"width_{suffix}")
    row["tss_status"] = tss.get("status")
    row["tss_enrichment"] = tss.get("enrichment")
    row["expected_motif_status"] = motif.get("status")
    row["best_motif_id"] = motif.get("best_motif_id")
    row["best_adjusted_p_value"] = motif.get("best_adjusted_p_value")
    row["ame_status"] = motif.get("ame_status")
    warnings = sample.get("warnings", [])
    row["warning_count"] = len(warnings) if isinstance(warnings, list) else 0
    return {key: row.get(key) for key in QC_SUMMARY_COLUMNS}


def _machine_tsv_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return format(value, ".15g")
    return str(value)


def _write_tsv(path: Path, columns: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: _machine_tsv_value(value)
                for key, value in row.items()
            })


def write_qc_summary_tsv(data: Mapping[str, object], path: Path) -> None:
    """Write one deterministic, flat QC summary row per derived sample."""
    _write_tsv(path, QC_SUMMARY_COLUMNS, [_summary_row(sample) for sample in _samples(data)])


def _project(mapping: Mapping[str, object], fields: Sequence[str]) -> dict[str, object]:
    return {field: mapping.get(field) for field in fields}


def _public_warning(
    warning: Mapping[str, object], *, default_sample_id: object = None
) -> dict[str, object]:
    return {
        "sample_id": warning.get("sample_id", default_sample_id),
        "family": warning.get("family", "general"),
        "status": warning.get("status", "missing"),
        "message": warning.get("message"),
    }


def _public_sample(sample: Mapping[str, object]) -> dict[str, object]:
    library = _nested(sample, "library")
    peak = _nested(sample, "peak")
    tss = _nested(sample, "tss")
    motif = _nested(sample, "motif")
    availability = sample.get("availability", {})
    public_availability = {
        family: {
            "status": (
                availability.get(family, {}).get("status")
                if isinstance(availability, Mapping)
                and isinstance(availability.get(family), Mapping)
                else "missing"
            ),
            "reason": (
                availability.get(family, {}).get("reason")
                if isinstance(availability, Mapping)
                and isinstance(availability.get(family), Mapping)
                else "not_evaluated"
            ),
        }
        for family in FAMILY_NAMES
    }
    insert_distribution = library.get("insert_size_distribution")
    public_library = _project(library, LIBRARY_FIELDS)
    public_library["insert_size_distribution"] = (
        [
            {
                "insert_size": row.get("insert_size"),
                "pair_count": row.get("pair_count"),
            }
            for row in insert_distribution
            if isinstance(row, Mapping)
        ]
        if isinstance(insert_distribution, list) else None
    )
    width_distribution = peak.get("width_distribution")
    fragments_per_peak_distribution = peak.get(
        "fragments_per_peak_distribution"
    )
    public_peak = _project(peak, PEAK_FIELDS)
    public_peak["width_distribution"] = (
        [
            {"width": row.get("width"), "peak_count": row.get("peak_count")}
            for row in width_distribution
            if isinstance(row, Mapping)
        ]
        if isinstance(width_distribution, list) else None
    )
    public_peak["fragments_per_peak_distribution"] = (
        [
            {
                "fragment_count": row.get("fragment_count"),
                "peak_count": row.get("peak_count"),
            }
            for row in fragments_per_peak_distribution
            if isinstance(row, Mapping)
        ]
        if isinstance(fragments_per_peak_distribution, list) else None
    )
    top_motifs = sample.get("top_motifs", [])
    warnings = sample.get("warnings", [])
    return {
        "sample_id": sample.get("sample_id"),
        "library_id": sample.get("library_id"),
        "input_group": sample.get("input_group"),
        "assay_target": sample.get("assay_target"),
        "is_control": sample.get("is_control"),
        "control_id": sample.get("control_id"),
        "expected_motif": sample.get("expected_motif"),
        "sample_kind": sample.get("sample_kind"),
        "demultiplex": _project(_nested(sample, "demultiplex"), DEMULTIPLEX_FIELDS),
        "library": public_library,
        "peak": public_peak,
        "tss": {
            "status": tss.get("status"),
            "enrichment": tss.get("enrichment"),
            "profile": tss.get("profile"),
            "score_status": tss.get("score_status"),
        },
        "motif": {
            "status": motif.get("status"),
            "best_motif_id": motif.get("best_motif_id"),
            "best_adjusted_p_value": motif.get("best_adjusted_p_value"),
            "ame_status": motif.get("ame_status"),
        },
        "availability": public_availability,
        "top_motifs": [
            _project(row, TOP_MOTIF_COLUMNS[3:])
            for row in top_motifs
            if isinstance(row, Mapping)
        ] if isinstance(top_motifs, list) else [],
        "warnings": [
            _public_warning(warning, default_sample_id=sample.get("sample_id"))
            for warning in warnings
            if isinstance(warning, Mapping)
        ] if isinstance(warnings, list) else [],
    }


def _json_payload(data: Mapping[str, object]) -> dict[str, object]:
    internal_samples = _samples(data)
    samples = [_public_sample(sample) for sample in internal_samples]
    raw_warnings = data.get("warnings", [])
    warnings = [
        _public_warning(warning)
        for warning in raw_warnings
        if isinstance(warning, Mapping)
    ] if isinstance(raw_warnings, list) else []
    target_count = sum(not bool(sample.get("is_control")) for sample in samples)
    return {
        "schema_version": data.get("schema_version", SCHEMA_VERSION),
        "generator_version": data.get("generator_version", GENERATOR_VERSION),
        "annotation_status": data.get("annotation_status"),
        "counts": {
            "samples": len(samples),
            "targets": target_count,
            "controls": len(samples) - target_count,
            "warnings": len(warnings),
        },
        "availability": _aggregate_availability(internal_samples),
        "metric_definitions": {
            "tss_enrichment": {
                "formula": "center_bin_signal / mean(terminal_100bp_flanks)",
                "before_bp": TSS_BEFORE_BP,
                "after_bp": TSS_BEFORE_BP,
                "bin_size_bp": TSS_BIN_SIZE,
                "center_position_bp": 0,
                "flank_bp_per_side": TSS_FLANK_BP,
            },
        },
        "samples": samples,
        "warnings": warnings,
    }


def write_qc_summary_json(data: Mapping[str, object], path: Path) -> None:
    """Write the structured report data with a stable, documented schema."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        json.dump(_json_payload(data), handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def _top_motif_rows(data: Mapping[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sample in _samples(data):
        if sample.get("is_control"):
            continue
        motifs = sample.get("top_motifs", [])
        if not isinstance(motifs, list):
            continue
        ordered = sorted(
            (motif for motif in motifs if isinstance(motif, Mapping)),
            key=lambda motif: (
                float(motif.get("rank", math.inf)) if motif.get("rank") is not None else math.inf,
                str(motif.get("motif_id", "")),
            ),
        )[:TOP_MOTIF_LIMIT]
        for motif in ordered:
            rows.append({
                "sample_id": sample.get("sample_id"),
                "assay_target": sample.get("assay_target"),
                "expected_motif": sample.get("expected_motif"),
                **{key: motif.get(key) for key in TOP_MOTIF_COLUMNS[3:]},
            })
    return rows


def write_top_motifs_tsv(data: Mapping[str, object], path: Path) -> None:
    """Write at most ten rank-sorted AME motifs for every target sample."""
    _write_tsv(path, TOP_MOTIF_COLUMNS, _top_motif_rows(data))


def _tss_profile_rows(data: Mapping[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sample in _samples(data):
        profile = _nested(sample, "tss").get("profile")
        if not isinstance(profile, Sequence) or isinstance(profile, (str, bytes)):
            continue
        parsed: list[tuple[int, object]] = []
        for item in profile:
            if not isinstance(item, Sequence) or len(item) != 2:
                continue
            parsed.append((int(item[0]), item[1]))
        for position_bp, signal in sorted(parsed):
            rows.append({
                "sample_id": sample.get("sample_id"),
                "position_bp": position_bp,
                "signal": signal,
            })
    return rows


def write_tss_profiles_tsv(data: Mapping[str, object], path: Path) -> None:
    """Write tidy relative-position TSS signals suitable for independent plotting."""
    _write_tsv(path, TSS_PROFILE_COLUMNS, _tss_profile_rows(data))


def render_table(
    columns, rows, *, empty_message, aria_label, row_limit: int | None = None
):
    """Render an escaped HTML table, including an explicit empty-state message."""
    all_rows = list(rows)
    if not all_rows:
        return f'<p class="empty">{html.escape(empty_message)}</p>'
    visible_rows = all_rows[:row_limit] if row_limit is not None else all_rows
    header = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(format_value(row.get(key)))}</td>"
            for key, _ in columns
        )
        + "</tr>"
        for row in visible_rows
    )
    table = (
        '<div class="table-scroll" role="region" '
        f'aria-label="{html.escape(aria_label)}" tabindex="0">'
        f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"
        "</div>"
    )
    if len(visible_rows) != len(all_rows):
        table += (
            '<p class="table-note">'
            f"Showing the first {len(visible_rows)} of {len(all_rows)} rows; "
            "the complete distribution is retained in qc_summary.json."
            "</p>"
        )
    return table


def render_bar_chart(
    title: str, rows: Sequence[Mapping[str, object]], *, value_key: str,
    label_key: str = "sample_id", axis_label: str = "Value",
) -> str:
    """Render a horizontally scrollable SVG whose bars retain readable labels."""
    numeric_rows: list[tuple[str, float, bool]] = []
    for row in rows:
        value = row.get(value_key)
        try:
            number = float(value) if value is not None else math.nan
        except (TypeError, ValueError):
            number = math.nan
        if math.isfinite(number):
            numeric_rows.append((str(row.get(label_key, "")), number, bool(row.get("is_control"))))
    if not numeric_rows:
        return f'<p class="empty">{html.escape("No numeric data available for " + title)}</p>'
    maximum = max((max(value, 0.0) for _, value, _ in numeric_rows), default=0.0) or 1.0
    longest_label = max(len(label) for label, _, _ in numeric_rows)
    slot_width = max(64.0, min(220.0, longest_label * 7.0 + 22.0))
    left, right, plot_bottom = 60.0, 24.0, 220.0
    height = max(330.0, plot_bottom + longest_label * 5.0 + 36.0)
    width = max(640.0, left + right + slot_width * len(numeric_rows))
    plot_height = plot_bottom - 32.0
    bar_width = slot_width * 0.64
    pattern_suffix = "".join(
        character.lower() if character.isalnum() else "-"
        for character in title
    ).strip("-") or "values"
    pattern_id = f"control-hatch-{pattern_suffix}"
    bars = []
    for index, (label, value, control) in enumerate(numeric_rows):
        x = left + index * slot_width + (slot_width - bar_width) / 2.0
        bar_height = max(0.0, value) / maximum * plot_height
        y = plot_bottom - bar_height
        fill = f"url(#{pattern_id})" if control else "#2563eb"
        stroke = "#334155" if control else "#2563eb"
        tooltip = f"{label}: {format_value(value)}"
        label_x = x + bar_width / 2.0
        label_y = plot_bottom + 18.0
        bars.append(
            f'<rect class="bar" x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" '
            f'height="{bar_height:.1f}" data-sample-kind="{"control" if control else "target"}" '
            f'fill="{fill}" stroke="{stroke}"><title>{html.escape(tooltip)}</title></rect>'
            f'<text x="{label_x:.1f}" y="{label_y:.1f}" '
            f'transform="rotate(45 {label_x:.1f} {label_y:.1f})" '
            f'text-anchor="start">{html.escape(label)}</text>'
        )
    return (
        '<div class="chart-scroll" role="region" '
        f'aria-label="{html.escape(title)} chart" tabindex="0">'
        f'<svg class="chart chart-wide" width="{width:.0f}" viewBox="0 0 {width:.0f} {height:.0f}" role="img" '
        f'aria-label="{html.escape(title)}"><title>{html.escape(title)}</title>'
        f'<defs><pattern id="{pattern_id}" width="6" height="6" patternUnits="userSpaceOnUse" '
        'patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="6" '
        'stroke="#334155" stroke-width="2"/></pattern></defs>'
        f'<text x="{left}" y="18" class="chart-title">{html.escape(title)}</text>'
        f'<text x="16" y="{plot_bottom / 2:.1f}" '
        f'transform="rotate(-90 16 {plot_bottom / 2:.1f})">{html.escape(axis_label)}</text>'
        f'<line x1="{left}" y1="{plot_bottom}" x2="{width - right}" y2="{plot_bottom}" class="axis"/>'
        + "".join(bars) + "</svg></div>"
    )


def _series_style(index: int, is_control: bool) -> tuple[str, str, str, str, str]:
    palette = ("#2563eb", "#7c3aed", "#0f766e", "#c2410c", "#be123c")
    return (
        f"S{index + 1:02d}",
        palette[index % len(palette)],
        "control" if is_control else "target",
        "square" if is_control else "circle",
        f"{index + 2} {index + 3}",
    )


def _series_marker(
    x: float, y: float, *, color: str, marker: str, tooltip: str
) -> str:
    if marker == "square":
        return (
            f'<rect x="{x - 3:.1f}" y="{y - 3:.1f}" width="6" height="6" '
            f'fill="{color}"><title>{html.escape(tooltip)}</title></rect>'
        )
    return (
        f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" '
        f'fill="{color}"><title>{html.escape(tooltip)}</title></circle>'
    )


def _series_legend_entry(
    series_key: str,
    sample_id: str,
    *,
    color: str,
    kind: str,
    marker: str,
    dash: str,
) -> str:
    swatch_marker = (
        f'<rect x="17" y="2" width="6" height="6" fill="{color}"/>'
        if marker == "square"
        else f'<circle cx="20" cy="5" r="3" fill="{color}"/>'
    )
    return (
        '<li>'
        '<svg class="series-swatch" viewBox="0 0 40 10" '
        'aria-hidden="true" focusable="false">'
        f'<line x1="0" y1="5" x2="40" y2="5" stroke="{color}" '
        f'stroke-width="2" stroke-dasharray="{dash}"/>{swatch_marker}</svg>'
        f'<span class="series-key">{series_key}</span>'
        f'<span class="series-label">{html.escape(sample_id)}</span>'
        f'<span class="series-kind"> ({kind}; {marker})</span>'
        '</li>'
    )


def _downsample_profile(
    profile: Sequence[tuple[float, float]], *, limit: int = 600
) -> list[tuple[float, float]]:
    """Retain local maxima when a series has more points than display pixels."""
    ordered = sorted(profile)
    if len(ordered) <= limit:
        return ordered
    chunk_size = math.ceil(len(ordered) / limit)
    sampled = [
        max(ordered[index:index + chunk_size], key=lambda point: point[1])
        for index in range(0, len(ordered), chunk_size)
    ]
    return sampled


def render_distribution_chart(
    title: str,
    rows: Sequence[Mapping[str, object]],
    *,
    x_key: str,
    value_key: str,
    x_axis_label: str,
    y_axis_label: str,
) -> str:
    """Render fixed-width numeric distributions with one trace per sample."""
    grouped: dict[str, tuple[bool, list[tuple[float, float]]]] = {}
    for row in rows:
        try:
            x_value = float(row.get(x_key))
            y_value = float(row.get(value_key))
        except (TypeError, ValueError):
            continue
        if not math.isfinite(x_value) or not math.isfinite(y_value):
            continue
        sample_id = str(row.get("sample_id", ""))
        if sample_id not in grouped:
            grouped[sample_id] = (bool(row.get("is_control")), [])
        grouped[sample_id][1].append((x_value, max(y_value, 0.0)))
    if not grouped:
        return (
            f'<p class="empty">{html.escape("No numeric data available for " + title)}</p>'
        )
    all_points = [
        point for _, profile in grouped.values() for point in profile
    ]
    x_min = min(point[0] for point in all_points)
    x_max = max(point[0] for point in all_points)
    y_max = max(point[1] for point in all_points) or 1.0
    x_span = x_max - x_min
    width, height, left, right, plot_bottom = 640, 280, 60, 24, 234
    plot_width = width - left - right
    plot_height = plot_bottom - 32
    def x_position(value: float) -> float:
        if x_span == 0:
            return left + plot_width / 2.0
        return left + (value - x_min) / x_span * plot_width

    axis_ticks = []
    x_tick_values = (
        (x_min,) if x_span == 0
        else (x_min, (x_min + x_max) / 2.0, x_max)
    )
    for value in x_tick_values:
        tick_x = x_position(value)
        axis_ticks.append(
            f'<line x1="{tick_x:.1f}" y1="{plot_bottom}" '
            f'x2="{tick_x:.1f}" y2="{plot_bottom + 5}" class="axis"/>'
            f'<text x="{tick_x:.1f}" y="{plot_bottom + 17}" '
            f'text-anchor="middle" class="axis-tick-label">'
            f'{html.escape(format_value(value))}</text>'
        )
    for value in (0.0, y_max / 2.0, y_max):
        y_position = plot_bottom - value / y_max * plot_height
        axis_ticks.append(
            f'<line x1="{left - 5}" y1="{y_position:.1f}" '
            f'x2="{left}" y2="{y_position:.1f}" class="axis"/>'
            f'<text x="{left - 8}" y="{y_position + 4:.1f}" '
            f'text-anchor="end" class="axis-tick-label">'
            f'{html.escape(format_value(value))}</text>'
        )
    traces = []
    legend = []
    for index, (sample_id, (is_control, raw_profile)) in enumerate(
        sorted(grouped.items())
    ):
        profile = _downsample_profile(raw_profile)
        plotted = [
            (
                x_position(x_value),
                plot_bottom - y_value / y_max * plot_height,
            )
            for x_value, y_value in profile
        ]
        coordinates = " ".join(f"{x:.1f},{y:.1f}" for x, y in plotted)
        series_key, color, kind, marker, dash = _series_style(index, is_control)
        tooltip = f"{sample_id}: {title} ({kind})"
        marker_x, marker_y = plotted[-1]
        traces.append(
            f'<polyline data-series-key="{series_key}" data-sample-kind="{kind}" '
            f'data-marker="{marker}" points="{coordinates}" fill="none" '
            f'stroke="{color}" stroke-width="2" stroke-dasharray="{dash}">'
            f'<title>{html.escape(tooltip)}</title></polyline>'
            + _series_marker(
                marker_x, marker_y, color=color, marker=marker, tooltip=tooltip
            )
        )
        legend.append(
            _series_legend_entry(
                series_key, sample_id, color=color, kind=kind,
                marker=marker, dash=dash,
            )
        )
    return (
        '<div class="line-chart">'
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="{html.escape(title)}"><title>{html.escape(title)}</title>'
        f'<text x="{left}" y="18" class="chart-title">{html.escape(title)}</text>'
        f'<text x="16" y="{plot_bottom / 2:.1f}" '
        f'transform="rotate(-90 16 {plot_bottom / 2:.1f})">'
        f'{html.escape(y_axis_label)}</text>'
        f'<text x="{width / 2:.1f}" y="{height - 6}" text-anchor="middle">'
        f'{html.escape(x_axis_label)}</text>'
        f'<line x1="{left}" y1="{plot_bottom}" x2="{width - right}" '
        f'y2="{plot_bottom}" class="axis"/>'
        + "".join(axis_ticks) + "".join(traces) + '</svg><ol class="series-legend">'
        + "".join(legend) + "</ol></div>"
    )


def render_line_chart(title: str, profiles: Mapping[str, object]) -> str:
    """Render a fixed plot with a unique textual and stroke key for each trace."""
    points: list[tuple[str, int, float, bool]] = []
    for sample_id, entry in profiles.items():
        if isinstance(entry, Mapping):
            profile = entry.get("profile", [])
            is_control = bool(entry.get("is_control"))
        else:
            profile = entry
            is_control = False
        if not isinstance(profile, Sequence) or isinstance(profile, (str, bytes)):
            continue
        for item in profile:
            if (
                not isinstance(item, Sequence)
                or isinstance(item, (str, bytes))
                or len(item) != 2
            ):
                continue
            position, signal = item
            try:
                numeric_signal = float(signal)
                numeric_position = int(position)
            except (TypeError, ValueError):
                continue
            if math.isfinite(numeric_signal):
                points.append((str(sample_id), numeric_position, numeric_signal, is_control))
    if not points:
        return f'<p class="empty">{html.escape("No TSS profile data available")}</p>'
    width, height, left, bottom = 640, 280, 60, 46
    x_min, x_max = min(point[1] for point in points), max(point[1] for point in points)
    y_max = max(max(point[2], 0.0) for point in points) or 1.0
    x_span = x_max - x_min or 1
    grouped: dict[str, tuple[bool, list[tuple[int, float]]]] = {}
    for sample_id, position, signal, is_control in points:
        if sample_id not in grouped:
            grouped[sample_id] = (is_control, [])
        grouped[sample_id][1].append((position, signal))
    plot_bottom = height - bottom
    plot_height = plot_bottom - 32
    paths = []
    legend = []
    for index, (sample_id, (is_control, profile)) in enumerate(sorted(grouped.items())):
        plotted = [
            (
                left + (position - x_min) / x_span * (width - left - 24),
                plot_bottom - max(signal, 0.0) / y_max * plot_height,
            )
            for position, signal in sorted(profile)
        ]
        coordinates = " ".join(f"{x:.1f},{y:.1f}" for x, y in plotted)
        series_key, color, kind, marker, dash = _series_style(index, is_control)
        tooltip = f"{sample_id}: TSS profile ({kind})"
        marker_x, marker_y = plotted[-1]
        paths.append(
            f'<polyline data-series-key="{series_key}" data-sample-kind="{kind}" '
            f'data-marker="{marker}" points="{coordinates}" fill="none" stroke="{color}" '
            f'stroke-width="2" stroke-dasharray="{dash}">'
            f'<title>{html.escape(tooltip)}</title></polyline>'
            + _series_marker(
                marker_x, marker_y, color=color, marker=marker, tooltip=tooltip
            )
        )
        legend.append(
            _series_legend_entry(
                series_key, sample_id, color=color, kind=kind,
                marker=marker, dash=dash,
            )
        )
    return (
        '<div class="line-chart">'
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(title)}">'
        f'<title>{html.escape(title)}</title><text x="{left}" y="18" class="chart-title">{html.escape(title)}</text>'
        f'<text x="16" y="{height / 2:.1f}" transform="rotate(-90 16 {height / 2:.1f})">Signal</text>'
        f'<text x="{width / 2:.1f}" y="{height - 6}" text-anchor="middle">Position relative to TSS (bp)</text>'
        f'<line x1="{left}" y1="{plot_bottom}" x2="{width - 18}" y2="{plot_bottom}" class="axis"/>'
        + "".join(paths) + '</svg><ol class="series-legend">'
        + "".join(legend) + "</ol></div>"
    )


def _section(section_id: str, title: str, content: str) -> str:
    return f'<section id="{section_id}"><h2>{html.escape(title)}</h2>{content}</section>'


def render_dashboard(data: Mapping[str, object]) -> str:
    """Return a self-contained descriptive dashboard with inline SVG and tables."""
    samples = _samples(data)
    summary_rows = [_summary_row(sample) for sample in samples]
    target_count = sum(not bool(sample.get("is_control")) for sample in samples)
    control_count = len(samples) - target_count
    counts = (
        '<p class="run-counts">'
        f'<strong>{len(samples)} samples</strong> · '
        f'<strong>{target_count} target{"s" if target_count != 1 else ""}</strong> · '
        f'<strong>{control_count} control{"s" if control_count != 1 else ""}</strong>'
        '</p>'
    )
    overview = render_table(
        [("sample_id", "Sample"), ("library_id", "Physical library"),
         ("input_group", "Input group"), ("assay_target", "Target"),
         ("sample_kind", "Sample kind"), ("control_id", "Matched control"),
         ("expected_motif", "Expected motif"), ("is_control", "IgG control"),
         ("annotation_status", "Annotation status")],
        [{**sample, "annotation_status": data.get("annotation_status")} for sample in samples],
        empty_message="No samples were supplied",
        aria_label="Run overview table",
    )
    aggregate_availability = _aggregate_availability(samples)
    availability_rows = [
        {
            "family": family.replace("_", " ").title(),
            **aggregate_availability[family],
        }
        for family in FAMILY_NAMES
    ]
    availability = render_table(
        [("family", "Input family"), ("status", "Run status"),
         ("computed", "Computed"), ("skipped", "Skipped"), ("empty", "Empty"),
         ("missing", "Missing"), ("failed", "Failed"),
         ("not_applicable", "Not applicable")],
        availability_rows,
        empty_message="No input-family availability was recorded",
        aria_label="Input-family availability table",
    )
    demux = render_bar_chart("Assigned read-pair fraction", summary_rows, value_key="assigned_fraction", axis_label="Fraction")
    demux += render_table(
        [("sample_id", "Sample"), ("total_read_pairs", "Total read pairs"),
         ("assigned_read_pairs", "Assigned read pairs"),
         ("ambiguous_read_pairs", "Ambiguous read pairs"),
         ("unassigned_read_pairs", "Unassigned read pairs"),
         ("assigned_fraction", "Assigned fraction"),
         ("ambiguous_fraction", "Ambiguous fraction"),
         ("unassigned_fraction", "Unassigned fraction"),
         ("sample_assigned_reads", "Sample assigned reads"),
         ("sample_assignment_fraction", "Sample assignment fraction")],
        summary_rows, empty_message="No demultiplexing metrics available",
        aria_label="Demultiplexing metrics table",
    )
    alignment = render_bar_chart("Mapped reads", summary_rows, value_key="mapped_percent", axis_label="Percent")
    alignment += render_table(
        [("sample_id", "Sample"), ("raw_total_reads", "Raw reads"),
         ("mapped_percent", "Mapped (%)"),
         ("properly_paired_percent", "Properly paired (%)"),
         ("mapq_filtered_reads", "MAPQ reads"),
         ("mapq_filtered_fragments", "MAPQ fragments"),
         ("mapq_filtered_fraction", "MAPQ fraction"),
         ("markdup_examined_reads", "Examined reads"),
         ("duplicate_total", "Duplicate reads"),
         ("duplicate_percent", "Duplicate (%)"),
         ("mitochondrial_percent", "Mitochondrial (%)"),
         ("estimated_library_size", "Estimated library size"),
         ("insert_size_total_pairs", "Insert pairs"),
         ("insert_size_min", "Insert min"), ("insert_size_q25", "Insert Q25"),
         ("insert_size_mean", "Insert mean"),
         ("insert_size_median", "Insert median"),
         ("insert_size_q75", "Insert Q75"), ("insert_size_max", "Insert max")],
        summary_rows, empty_message="No alignment metrics available",
        aria_label="Alignment and library QC table",
    )
    insert_rows: list[dict[str, object]] = []
    for sample in samples:
        distribution = _nested(sample, "library").get("insert_size_distribution")
        if not isinstance(distribution, list):
            continue
        for row in distribution:
            if isinstance(row, Mapping):
                insert_rows.append({
                    "sample_id": sample.get("sample_id"),
                    "is_control": sample.get("is_control"),
                    "insert_size": row.get("insert_size"),
                    "pair_count": row.get("pair_count"),
                })
    insert = render_distribution_chart(
        "Insert-size distribution", insert_rows,
        x_key="insert_size", value_key="pair_count",
        x_axis_label="Insert size (bp)", y_axis_label="Read pairs",
    )
    insert += render_table(
        [("sample_id", "Sample"), ("insert_size", "Insert size (bp)"),
         ("pair_count", "Read pairs")],
        insert_rows, empty_message="No insert-size distribution data available",
        aria_label="Insert-size distribution table",
        row_limit=2000,
    )
    peak_rows = [row for row in summary_rows if not row.get("is_control")]
    peaks = render_bar_chart("FRiP", peak_rows, value_key="frip", axis_label="FRiP")
    peaks += render_table(
        [("sample_id", "Sample"), ("peak_count", "Peak count"),
         ("total_covered_bases", "Covered bases"),
         ("peak_width_min", "Width min"), ("peak_width_q25", "Width Q25"),
         ("peak_width_mean", "Width mean"),
         ("peak_width_median", "Width median"),
         ("peak_width_q75", "Width Q75"), ("peak_width_max", "Width max"),
         ("total_fragments", "Usable fragments"),
         ("fragments_in_peaks", "Fragments in peaks"), ("frip", "FRiP")],
        summary_rows, empty_message="No peak metrics available",
        aria_label="Peaks and FRiP table",
    )
    peak_width_rows: list[dict[str, object]] = []
    for sample in samples:
        distribution = _nested(sample, "peak").get("width_distribution")
        if not isinstance(distribution, list):
            continue
        for row in distribution:
            if isinstance(row, Mapping):
                peak_width_rows.append({
                    "sample_id": sample.get("sample_id"),
                    "is_control": sample.get("is_control"),
                    "width": row.get("width"),
                    "peak_count": row.get("peak_count"),
                })
    peak_width = render_distribution_chart(
        "Peak-width distribution", peak_width_rows,
        x_key="width", value_key="peak_count",
        x_axis_label="Peak width (bp)", y_axis_label="Peaks",
    )
    peak_width += render_table(
        [("sample_id", "Sample"), ("width", "Peak width (bp)"),
         ("peak_count", "Peaks")],
        peak_width_rows, empty_message="No peak-width distribution data available",
        aria_label="Peak-width distribution table",
        row_limit=2000,
    )
    profiles = {
        str(sample.get("sample_id", "")): {
            "profile": _nested(sample, "tss").get("profile", []),
            "is_control": sample.get("is_control"),
        }
        for sample in samples if isinstance(_nested(sample, "tss").get("profile"), Sequence)
    }
    tss = render_line_chart("TSS profiles", profiles)
    tss += render_table(
        [("sample_id", "Sample"), ("tss_status", "Status"), ("tss_enrichment", "TSS enrichment")],
        summary_rows, empty_message="No TSS metrics available",
        aria_label="TSS enrichment table",
    )
    motif_rows = _top_motif_rows(data)
    expected = render_table(
        [("sample_id", "Sample"), ("expected_motif", "Expected motif"),
         ("expected_motif_status", "Expected motif status"),
         ("best_motif_id", "Best motif"),
         ("best_adjusted_p_value", "Best adjusted significance"),
         ("ame_status", "AME status")],
        peak_rows, empty_message="No target motif metrics available",
        aria_label="Expected motif table",
    )
    expected += "<h3>Top AME motifs</h3>" + render_table(
        [("sample_id", "Sample"), ("rank", "Rank"), ("motif_id", "Motif"),
         ("motif_alt_id", "Alternate ID"), ("adjusted_p_value", "Adjusted p-value")],
        motif_rows, empty_message="No AME motifs available",
        aria_label="Top AME motifs table",
    )
    warnings = data.get("warnings", [])
    warning_rows = warnings if isinstance(warnings, list) else []
    warning_section = render_table(
        [("sample_id", "Sample"), ("family", "Input family"),
         ("status", "Status"), ("message", "Warning")],
        warning_rows, empty_message="No warnings recorded",
        aria_label="Warnings table",
    )
    body = "".join((
        _section("run-overview", "Run overview", counts + overview),
        _section("input-availability", "Input-family availability", availability),
        '<p class="legend">Bar charts — Targets (solid); IgG controls (outlined) '
        'and hatched. Line charts use the per-series swatches.</p>',
        _section("demultiplexing", "Demultiplexing", demux),
        _section("alignment", "Alignment and library QC", alignment),
        _section("insert-size-distribution", "Insert-size distribution", insert),
        _section("peaks-frip", "Peaks and FRiP", peaks),
        _section("peak-width-distribution", "Peak-width distribution", peak_width),
        _section("tss-enrichment", "TSS enrichment", tss),
        _section("motif-enrichment", "Motif enrichment", expected),
        _section("warnings", "Warnings", warning_section),
    ))
    return """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Consolidated QC dashboard</title><style>
body{font-family:system-ui,sans-serif;line-height:1.45;margin:0;color:#172033;background:#f8fafc}main{max-width:1100px;margin:auto;padding:1.5rem}section{background:#fff;border:1px solid #dbe3ee;border-radius:.5rem;padding:1rem;margin:1rem 0}h1,h2,h3{margin-top:0}.legend,.empty{color:#475569}.table-scroll{max-width:100%;overflow-x:auto}table{border-collapse:collapse;width:100%;margin:.75rem 0}th,td{border:1px solid #dbe3ee;padding:.35rem;text-align:left;vertical-align:top}th{background:#eff6ff}.chart{width:100%;height:auto;background:#fff}.axis{stroke:#64748b}.chart-title{font-weight:700}
.chart-scroll{max-width:100%;overflow-x:auto}.chart-wide{width:auto;min-width:100%;max-width:none}.series-legend{display:grid;grid-template-columns:repeat(auto-fit,minmax(18rem,1fr));gap:.25rem 1rem;padding-left:1.5rem}.series-swatch{width:2.5rem;height:.75rem;vertical-align:middle;margin-right:.35rem}.series-key{display:inline-block;min-width:2.5rem;font-weight:700}.series-label{font-family:ui-monospace,monospace}.run-counts{font-size:1.05rem}.table-note{color:#475569;font-size:.9rem}
</style></head><body><main><h1>Consolidated QC dashboard</h1><p>Descriptive technical and biological QC summary; no biological thresholds are applied.</p>""" + body + "</main></body></html>"


def write_outputs(data: Mapping[str, object], outdir: Path) -> None:
    """Write all dashboard artifacts to an existing or newly created directory."""
    outdir.mkdir(parents=True, exist_ok=True)
    write_qc_summary_tsv(data, outdir / "qc_summary.tsv")
    write_qc_summary_json(data, outdir / "qc_summary.json")
    write_top_motifs_tsv(data, outdir / "top_motifs.tsv")
    write_tss_profiles_tsv(data, outdir / "tss_profiles.tsv")
    (outdir / "qc_dashboard.html").write_text(render_dashboard(data), encoding="utf-8")


def _validate_output_files(outdir: Path) -> None:
    for name in DASHBOARD_OUTPUT_FILENAMES:
        if not (outdir / name).is_file() or (outdir / name).stat().st_size == 0:
            raise DashboardInputError(f"dashboard output {name} is empty")


def _replace_path(source: Path, destination: Path) -> None:
    """Replace one complete filesystem path; kept as a fault-injection seam."""
    os.replace(source, destination)


def _remove_generated_path(path: Path) -> None:
    """Remove a generated file or directory without following directory symlinks."""
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


def _validate_publication_target(outdir: Path) -> None:
    """Reject directory swaps that could remove files outside this bundle."""
    resolved = outdir.resolve()
    if resolved == Path.cwd().resolve() or resolved == Path(resolved.anchor):
        raise DashboardInputError(
            "dashboard output must be a dedicated dashboard output directory"
        )
    if outdir.is_symlink() or (outdir.exists() and not outdir.is_dir()):
        raise DashboardInputError(
            "dashboard output must be a dedicated dashboard output directory"
        )
    if not outdir.exists():
        return
    allowed = set(DASHBOARD_OUTPUT_FILENAMES)
    unexpected = sorted(
        child.name
        for child in outdir.iterdir()
        if child.name not in allowed or not child.is_file() or child.is_symlink()
    )
    if unexpected:
        raise DashboardInputError(
            "dashboard output must be a dedicated dashboard output directory; "
            f"unexpected entry {unexpected[0]}"
        )


def write_outputs_atomically(data: Mapping[str, object], outdir: Path) -> None:
    """Publish the validated five-file bundle as one rollback-safe directory."""
    _validate_publication_target(outdir)
    outdir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{outdir.name}.tmp-", dir=outdir.parent))
    backup = Path(tempfile.mkdtemp(prefix=f".{outdir.name}.backup-", dir=outdir.parent))
    backup.rmdir()
    had_existing_output = outdir.exists() or outdir.is_symlink()
    try:
        write_outputs(data, temporary)
        _validate_output_files(temporary)
        if had_existing_output:
            _replace_path(outdir, backup)
        _replace_path(temporary, outdir)
    except BaseException as publish_error:
        try:
            if backup.exists() or backup.is_symlink():
                if outdir.exists() or outdir.is_symlink():
                    _remove_generated_path(outdir)
                _replace_path(backup, outdir)
            elif not had_existing_output and (outdir.exists() or outdir.is_symlink()):
                _remove_generated_path(outdir)
        except BaseException as rollback_error:
            publish_error.add_note(
                f"dashboard rollback also failed; prior output is retained at {backup}: "
                f"{rollback_error}"
            )
        raise
    else:
        if backup.exists() or backup.is_symlink():
            _remove_generated_path(backup)
    finally:
        if temporary.exists() or temporary.is_symlink():
            _remove_generated_path(temporary)


def _input_paths(directory: Path, pattern: str, *, label: str) -> list[Path]:
    if not directory.is_dir():
        raise DashboardInputError(f"{label} directory does not exist: {directory}")
    return sorted(path for path in directory.rglob(pattern) if path.is_file())


def _read_status_rows(paths: Sequence[Path], *, label: str) -> dict[str, dict[str, object]]:
    """Read one-row status TSVs while retaining textual status values."""
    records: dict[str, dict[str, object]] = {}
    for path in paths:
        headers, rows = _read_tsv(path, label=label)
        if "sample_id" not in headers or len(rows) != 1:
            raise DashboardInputError(f"{path}: {label} must contain exactly one sample_id row")
        raw = rows[0]
        sample_id = _required_text(raw, "sample_id", label=label)
        if sample_id in records:
            raise DashboardInputError(f"duplicate sample_id {sample_id}")
        record: dict[str, object] = {"sample_id": sample_id}
        for key, value in raw.items():
            if key == "sample_id":
                continue
            record[key] = None if not value.strip() else value.strip()
        records[sample_id] = record
    return dict(sorted(records.items()))


def _read_motif_metrics(paths: Sequence[Path]) -> dict[str, dict[str, object]]:
    records = _read_status_rows(paths, label="motif metrics")
    for record in records.values():
        status = record.pop("expected_motif_status", record.get("status"))
        if "status" in record:
            record.pop("status")
        record["status"] = status
        value = record.get("best_adjusted_p_value")
        if value is not None:
            record["best_adjusted_p_value"] = finite_number(
                value, label="best_adjusted_p_value", minimum=0
            )
    return records


def _read_tss_metrics(directory: Path) -> dict[str, dict[str, object]]:
    statuses = _read_status_rows(
        _input_paths(directory, "*.tss_status.tsv", label="TSS"), label="TSS status"
    )
    profiles: dict[str, list[tuple[int, float | None]]] = {}
    for path in _input_paths(directory, "*.tss_profile.tsv", label="TSS"):
        suffix = ".tss_profile.tsv"
        sample_id = path.name[:-len(suffix)]
        if not sample_id or sample_id in profiles:
            raise DashboardInputError(f"duplicate TSS profile sample_id {sample_id}")
        profiles[sample_id] = read_tss_profile(path)
    records: dict[str, dict[str, object]] = {}
    for sample_id in sorted(set(statuses) | set(profiles)):
        status_record = statuses.get(sample_id)
        profile = profiles.get(sample_id)
        if status_record is None:
            score, score_status = _tss_enrichment_result(profile or [])
            records[sample_id] = {
                "sample_id": sample_id,
                "status": "failed",
                "reason": "status_missing",
                "profile": profile,
                "enrichment": score,
                "score_status": score_status,
            }
            continue
        record = dict(status_record)
        if profile is None:
            if record.get("status") == "computed":
                record["status"] = "failed"
                record["reason"] = "profile_missing"
            record["profile"] = None
            record["enrichment"] = None
            record["score_status"] = "profile_missing"
        else:
            score, score_status = _tss_enrichment_result(profile)
            record["profile"] = profile
            record["enrichment"] = score
            record["score_status"] = score_status
        records[sample_id] = record
    return dict(sorted(records.items()))


def _sample_id_for_ame(path: Path) -> str:
    """Recover the sample directory name from the AME producer's ``sample/ame`` tree."""
    return path.parent.parent.name if path.parent.name == "ame" else path.parent.name


def _read_ame_motifs_from_directory(
    directory: Path,
) -> tuple[
    dict[str, list[dict[str, object]]],
    dict[str, list[dict[str, object]]],
]:
    """Read each AME result once and return complete and bounded mappings."""
    complete: dict[str, list[dict[str, object]]] = {}
    for path in _input_paths(directory, "ame.tsv", label="AME"):
        sample_id = _sample_id_for_ame(path)
        if not sample_id or sample_id in complete:
            raise DashboardInputError(f"duplicate AME sample_id {sample_id}")
        complete[sample_id] = read_all_ame(path)
    complete = dict(sorted(complete.items()))
    top = {
        sample_id: records[:TOP_MOTIF_LIMIT]
        for sample_id, records in complete.items()
    }
    return complete, top


def _read_top_motifs_from_directory(directory: Path) -> dict[str, list[dict[str, object]]]:
    """Compatibility wrapper returning the bounded AME mapping."""
    return _read_ame_motifs_from_directory(directory)[1]


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="qc_dashboard.py")
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--demux-dir", required=True, type=Path)
    parser.add_argument("--library-dir", required=True, type=Path)
    parser.add_argument("--insert-dir", required=True, type=Path)
    parser.add_argument("--peak-dir", required=True, type=Path)
    parser.add_argument("--tss-dir", required=True, type=Path)
    parser.add_argument("--motif-dir", required=True, type=Path)
    parser.add_argument("--ame-dir", required=True, type=Path)
    parser.add_argument(
        "--annotation-status", required=True,
        choices=("skipped_no_annotation", "computed_bed", "computed_gtf"),
    )
    parser.add_argument(
        "--motif-analysis-status",
        choices=("auto", "skipped_no_database", "computed"),
        default="auto",
    )
    parser.add_argument("--outdir", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Build all five dashboard files, returning a shell-compatible status code."""
    args = _build_argument_parser().parse_args(argv)
    try:
        metadata = load_metadata(args.metadata)
        demultiplex = read_demultiplex_metrics(
            _input_paths(args.demux_dir, "*.json", label="demultiplex"), metadata
        )
        libraries = read_library_metrics(
            _input_paths(args.library_dir, "*.library_qc.tsv", label="library")
        )
        insert_sizes = read_insert_size_distributions(
            _input_paths(
                args.insert_dir, "*.insert_size_distribution.tsv", label="insert"
            )
        )
        peaks = read_peak_metrics(_input_paths(args.peak_dir, "*.peak_qc.tsv", label="peak"))
        peak_widths = read_peak_width_distributions(
            _input_paths(
                args.peak_dir, "*.peak_qc.width_histogram.tsv", label="peak"
            )
        )
        fragments_per_peak = read_fragments_per_peak_distributions(
            _input_paths(
                args.peak_dir,
                "*.peak_qc.fragments_per_peak.tsv",
                label="peak",
            )
        )
        tss = _read_tss_metrics(args.tss_dir)
        motifs = _read_motif_metrics(
            _input_paths(args.motif_dir, "*.motif_qc.tsv", label="motif")
        )
        ame_motifs, top_motifs = _read_ame_motifs_from_directory(args.ame_dir)
        data = build_report_data(
            metadata, demultiplex, libraries, peaks, tss, motifs, top_motifs,
            annotation_status=args.annotation_status,
            insert_sizes=insert_sizes,
            peak_widths=peak_widths,
            fragments_per_peak=fragments_per_peak,
            ame_motifs=ame_motifs,
            motif_analysis_status=args.motif_analysis_status,
        )
        write_outputs_atomically(data, args.outdir)
    except DashboardInputError as error:
        print(f"qc_dashboard.py: error: {error}", file=__import__("sys").stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
