#!/usr/bin/env python3
"""Parse and join the QC inputs used by the consolidated dashboard."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import statistics
import tempfile
from pathlib import Path
from typing import Mapping, Sequence

from motif_qc import read_ame


SCHEMA_VERSION = 1
TSS_BEFORE_BP = 3000
TSS_BIN_SIZE = 10
TSS_FLANK_BP = 100
TOP_MOTIF_LIMIT = 10

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
    "fragments_in_peaks", "frip", "tss_status", "tss_enrichment",
    "expected_motif_status", "best_motif_id", "best_adjusted_p_value",
    "ame_status", "warning_count",
]

TOP_MOTIF_COLUMNS = [
    "sample_id", "assay_target", "expected_motif", "rank", "motif_id",
    "motif_alt_id", "adjusted_p_value", "p_value", "effect",
    "positive_sequences",
]
TSS_PROFILE_COLUMNS = ["sample_id", "position_bp", "signal"]


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


def _read_single_sample_tsv(paths: Sequence[Path], *, label: str) -> dict[str, dict[str, object]]:
    parsed: dict[str, dict[str, object]] = {}
    for path in paths:
        headers, rows = _read_tsv(path, label=label)
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
            record[name] = None if not value.strip() else finite_number(
                value, label=f"{sample_id} {name}", minimum=0
            )
        parsed[sample_id] = record
    return dict(sorted(parsed.items()))


def read_library_metrics(paths: Sequence[Path]) -> dict[str, dict[str, object]]:
    """Read one library-QC summary row for every derived sample."""
    return _read_single_sample_tsv(paths, label="library metrics")


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
        sample_id = raw_metrics.get("sample_id", "").strip()
        if not sample_id:
            raise DashboardInputError(f"{path}: peak sample_id is required")
        if sample_id in parsed:
            raise DashboardInputError(f"duplicate sample_id {sample_id}")
        record: dict[str, object] = {"sample_id": sample_id}
        for metric, value in raw_metrics.items():
            if metric != "sample_id":
                record[metric] = None if not value.strip() else finite_number(
                    value, label=f"{sample_id} peak {metric}", minimum=0
                )
        parsed[sample_id] = record
    return dict(sorted(parsed.items()))


def read_tss_profile(
    path: Path, *, before_bp: int = TSS_BEFORE_BP, bin_size: int = TSS_BIN_SIZE
) -> list[tuple[int, float]]:
    """Read the one aggregate deepTools profile row with a fixed bin count."""
    if before_bp <= 0 or bin_size <= 0 or (2 * before_bp) % bin_size:
        raise DashboardInputError("TSS profile window must divide evenly into bins")
    expected_bins = (2 * before_bp) // bin_size
    candidates: list[list[float]] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise DashboardInputError(f"cannot read TSS profile {path}: {error}") from error
    for line_number, raw_line in enumerate(lines, start=1):
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        values: list[float] = []
        for field in reversed(raw_line.split("\t")):
            try:
                value = float(field)
            except ValueError:
                break
            if not math.isfinite(value):
                raise DashboardInputError(f"{path}:{line_number}: TSS profile values must be finite")
            values.append(value)
        values.reverse()
        if len(values) == expected_bins:
            candidates.append(values)
    if not candidates:
        raise DashboardInputError(f"{path}: no TSS profile row with {expected_bins} numeric bins")
    if len(candidates) != 1:
        raise DashboardInputError(f"{path}: TSS profile has multiple compatible data rows")
    return [(-before_bp + index * bin_size, value) for index, value in enumerate(candidates[0])]


def calculate_tss_enrichment(
    profile: Sequence[tuple[int, float]], *, flank_bp: int = TSS_FLANK_BP
) -> float | None:
    """Return center signal divided by the two terminal TSS flank means."""
    if flank_bp <= 0 or flank_bp % TSS_BIN_SIZE:
        raise DashboardInputError("TSS flank_bp must be a positive multiple of bin size")
    flank_bins = flank_bp // TSS_BIN_SIZE
    if len(profile) < 2 * flank_bins + 1:
        raise DashboardInputError("TSS profile is too short for the requested flanks")
    center_value = finite_number(profile[len(profile) // 2][1], label="TSS center value")
    flank_values = [finite_number(value, label="TSS flank value") for _, value in profile[:flank_bins]]
    flank_values += [finite_number(value, label="TSS flank value") for _, value in profile[-flank_bins:]]
    flank_mean = statistics.fmean(flank_values)
    return None if flank_mean == 0 else center_value / flank_mean


def read_top_ame(path: Path, *, limit: int = TOP_MOTIF_LIMIT) -> list[dict[str, object]]:
    """Return up to ``limit`` AME records in deterministic significance order."""
    if limit < 0:
        raise DashboardInputError("AME motif limit must be non-negative")
    try:
        records = read_ame(path)
    except ValueError as error:
        raise DashboardInputError(str(error)) from error
    top = [record for record in records if record.motif_id != "__NO_PEAKS__"]
    top.sort(key=lambda record: (
        record.adjusted_p_value,
        record.rank if record.rank is not None else math.inf,
        record.motif_id,
    ))
    return [
        {
            "motif_id": record.motif_id,
            "motif_alt_id": record.motif_alt_id,
            "adjusted_p_value": record.adjusted_p_value,
            "p_value": record.p_value,
            "effect": record.effect,
            "positive_sequences": record.positive_sequences,
            "rank": record.rank,
        }
        for record in top[:limit]
    ]


def _empty_sample(record: Mapping[str, object]) -> dict[str, object]:
    return {
        **dict(record),
        "sample_kind": "control" if record["is_control"] else "target",
        "demultiplex": {
            "total_read_pairs": None, "assigned_read_pairs": None,
            "ambiguous_read_pairs": None, "unassigned_read_pairs": None,
            "assigned_fraction": None, "ambiguous_fraction": None,
            "unassigned_fraction": None, "sample_assigned_reads": None,
            "sample_assignment_fraction": None,
        },
        "library": {},
        "peak": {"frip": None},
        "tss": {"status": None, "enrichment": None, "profile": None},
        "motif": {"status": "not_applicable_control" if record["is_control"] else None},
        "top_motifs": [],
        "warnings": [],
    }


def build_report_data(
    metadata: Mapping[str, Mapping[str, object]],
    demultiplex: Mapping[str, Mapping[str, object]],
    libraries: Mapping[str, Mapping[str, object]],
    peaks: Mapping[str, Mapping[str, object]],
    tss: Mapping[str, Mapping[str, object]],
    motifs: Mapping[str, Mapping[str, object]],
    top_motifs: Mapping[str, Sequence[Mapping[str, object]]],
    annotation_status: str,
) -> dict[str, object]:
    """Initialize optional QC families, then overlay validated sample metrics."""
    samples_by_id = {sample_id: _empty_sample(record) for sample_id, record in metadata.items()}
    sample_ids = set(samples_by_id)
    for family_name, family in (("library", libraries), ("peak", peaks), ("tss", tss),
                                ("motif", motifs), ("top motifs", top_motifs)):
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
            sample["warnings"].append({"message": f"{sample_id}: missing demultiplex metrics for {library_id}"})
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
        if sample_id in libraries:
            sample["library"].update(libraries[sample_id])
        else:
            sample["warnings"].append({"message": f"{sample_id}: missing library metrics"})
        if sample_id in peaks:
            sample["peak"].update(peaks[sample_id])
        elif not sample["is_control"]:
            sample["warnings"].append({"message": f"{sample_id}: missing peak metrics"})
        if sample_id in tss:
            sample["tss"].update(tss[sample_id])
        else:
            sample["warnings"].append({"message": f"{sample_id}: missing TSS metrics"})
        if not sample["is_control"]:
            if sample_id in motifs:
                sample["motif"].update(motifs[sample_id])
            else:
                sample["warnings"].append({"message": f"{sample_id}: missing motif metrics"})
            sample["top_motifs"] = [dict(item) for item in top_motifs.get(sample_id, ())]
    samples = [samples_by_id[sample_id] for sample_id in sorted(samples_by_id)]
    warnings = [warning for sample in samples for warning in sample["warnings"]]
    return {
        "schema_version": SCHEMA_VERSION,
        "annotation_status": annotation_status,
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
    row["tss_status"] = tss.get("status")
    row["tss_enrichment"] = tss.get("enrichment")
    row["expected_motif_status"] = motif.get("status")
    row["best_motif_id"] = motif.get("best_motif_id")
    row["best_adjusted_p_value"] = motif.get("best_adjusted_p_value")
    row["ame_status"] = motif.get("ame_status")
    warnings = sample.get("warnings", [])
    row["warning_count"] = len(warnings) if isinstance(warnings, list) else 0
    return {key: row.get(key) for key in QC_SUMMARY_COLUMNS}


def _write_tsv(path: Path, columns: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if value is None else value for key, value in row.items()})


def write_qc_summary_tsv(data: Mapping[str, object], path: Path) -> None:
    """Write one deterministic, flat QC summary row per derived sample."""
    _write_tsv(path, QC_SUMMARY_COLUMNS, [_summary_row(sample) for sample in _samples(data)])


def _json_payload(data: Mapping[str, object]) -> dict[str, object]:
    samples = [dict(sample) for sample in _samples(data)]
    warnings = data.get("warnings", [])
    return {
        "schema_version": data.get("schema_version", SCHEMA_VERSION),
        "annotation_status": data.get("annotation_status"),
        "sample_count": len(samples),
        "warning_count": len(warnings) if isinstance(warnings, list) else 0,
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
                float(motif.get("adjusted_p_value", math.inf))
                if motif.get("adjusted_p_value") is not None else math.inf,
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


def render_table(columns, rows, *, empty_message):
    """Render an escaped HTML table, including an explicit empty-state message."""
    if not rows:
        return f'<p class="empty">{html.escape(empty_message)}</p>'
    header = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(format_value(row.get(key)))}</td>"
            for key, _ in columns
        )
        + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"


def render_bar_chart(
    title: str, rows: Sequence[Mapping[str, object]], *, value_key: str,
    label_key: str = "sample_id", axis_label: str = "Value",
) -> str:
    """Render a compact, accessible inline SVG bar chart without remote assets."""
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
    maximum = max(value for _, value, _ in numeric_rows) or 1.0
    width, height, left, bottom = 640, 260, 60, 46
    plot_height = height - bottom - 32
    bar_width = max(12, (width - left - 24) / len(numeric_rows) * 0.7)
    bars = []
    for index, (label, value, control) in enumerate(numeric_rows):
        x = left + (index + 0.5) * (width - left - 24) / len(numeric_rows) - bar_width / 2
        bar_height = max(0.0, value / maximum * plot_height)
        y = height - bottom - bar_height
        fill = "url(#control-hatch)" if control else "#2563eb"
        stroke = "#334155" if control else "#2563eb"
        tooltip = f"{label}: {format_value(value)}"
        bars.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" '
            f'fill="{fill}" stroke="{stroke}"><title>{html.escape(tooltip)}</title></rect>'
            f'<text x="{x + bar_width / 2:.1f}" y="{height - 24}" text-anchor="middle">{html.escape(label)}</text>'
        )
    return (
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="{html.escape(title)}"><title>{html.escape(title)}</title>'
        '<defs><pattern id="control-hatch" width="6" height="6" patternUnits="userSpaceOnUse" '
        'patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="6" '
        'stroke="#334155" stroke-width="2"/></pattern></defs>'
        f'<text x="{left}" y="18" class="chart-title">{html.escape(title)}</text>'
        f'<text x="16" y="{height / 2:.1f}" transform="rotate(-90 16 {height / 2:.1f})">{html.escape(axis_label)}</text>'
        f'<line x1="{left}" y1="{height - bottom}" x2="{width - 18}" y2="{height - bottom}" class="axis"/>'
        + "".join(bars) + "</svg>"
    )


def render_line_chart(title: str, profiles: Mapping[str, Sequence[tuple[int, object]]]) -> str:
    """Render escaped profile labels and SVG title tooltips for TSS traces."""
    points = [
        (str(sample_id), int(position), float(signal))
        for sample_id, profile in profiles.items()
        for position, signal in profile
        if math.isfinite(float(signal))
    ]
    if not points:
        return f'<p class="empty">{html.escape("No TSS profile data available")}</p>'
    width, height, left, bottom = 640, 260, 60, 46
    x_min, x_max = min(point[1] for point in points), max(point[1] for point in points)
    y_max = max(point[2] for point in points) or 1.0
    x_span = x_max - x_min or 1
    grouped: dict[str, list[tuple[int, float]]] = {}
    for sample_id, position, signal in points:
        grouped.setdefault(sample_id, []).append((position, signal))
    palette = ("#2563eb", "#7c3aed", "#0f766e", "#c2410c", "#be123c")
    paths = []
    for index, (sample_id, profile) in enumerate(sorted(grouped.items())):
        coordinates = " ".join(
            f"{left + (position - x_min) / x_span * (width - left - 24):.1f},{height - bottom - signal / y_max * (height - bottom - 32):.1f}"
            for position, signal in sorted(profile)
        )
        color = palette[index % len(palette)]
        tooltip = f"{sample_id}: TSS profile"
        paths.append(
            f'<polyline points="{coordinates}" fill="none" stroke="{color}" stroke-width="2">'
            f'<title>{html.escape(tooltip)}</title></polyline>'
            f'<text x="{width - 18}" y="{36 + index * 16}" text-anchor="end" fill="{color}">{html.escape(sample_id)}</text>'
        )
    return (
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(title)}">'
        f'<title>{html.escape(title)}</title><text x="{left}" y="18" class="chart-title">{html.escape(title)}</text>'
        f'<text x="16" y="{height / 2:.1f}" transform="rotate(-90 16 {height / 2:.1f})">Signal</text>'
        f'<text x="{width / 2:.1f}" y="{height - 6}" text-anchor="middle">Position relative to TSS (bp)</text>'
        f'<line x1="{left}" y1="{height - bottom}" x2="{width - 18}" y2="{height - bottom}" class="axis"/>'
        + "".join(paths) + "</svg>"
    )


def _section(section_id: str, title: str, content: str) -> str:
    return f'<section id="{section_id}"><h2>{html.escape(title)}</h2>{content}</section>'


def render_dashboard(data: Mapping[str, object]) -> str:
    """Return a self-contained descriptive dashboard with inline SVG and tables."""
    samples = _samples(data)
    summary_rows = [_summary_row(sample) for sample in samples]
    overview = render_table(
        [("sample_id", "Sample"), ("assay_target", "Target"), ("is_control", "IgG control"),
         ("annotation_status", "Annotation status")],
        [{**sample, "annotation_status": data.get("annotation_status")} for sample in samples],
        empty_message="No samples were supplied",
    )
    demux = render_bar_chart("Assigned read-pair fraction", summary_rows, value_key="assigned_fraction", axis_label="Fraction")
    demux += render_table(
        [("sample_id", "Sample"), ("assigned_fraction", "Assigned fraction"),
         ("ambiguous_fraction", "Ambiguous fraction"), ("unassigned_fraction", "Unassigned fraction")],
        summary_rows, empty_message="No demultiplexing metrics available",
    )
    alignment = render_bar_chart("Mapped reads", summary_rows, value_key="mapped_percent", axis_label="Percent")
    alignment += render_table(
        [("sample_id", "Sample"), ("mapped_percent", "Mapped (%)"),
         ("properly_paired_percent", "Properly paired (%)")],
        summary_rows, empty_message="No alignment metrics available",
    )
    peak_rows = [row for row in summary_rows if not row.get("is_control")]
    peaks = render_bar_chart("FRiP", peak_rows, value_key="frip", axis_label="FRiP")
    peaks += render_table(
        [("sample_id", "Sample"), ("peak_count", "Peak count"), ("frip", "FRiP")],
        peak_rows, empty_message="No target peak metrics available",
    )
    profiles = {
        str(sample.get("sample_id", "")): _nested(sample, "tss").get("profile", [])
        for sample in samples if isinstance(_nested(sample, "tss").get("profile"), Sequence)
    }
    safe_profiles = {
        sample_id: profile for sample_id, profile in profiles.items()
        if not isinstance(profile, (str, bytes))
    }
    tss = render_line_chart("TSS profiles", safe_profiles)
    tss += render_table(
        [("sample_id", "Sample"), ("tss_status", "Status"), ("tss_enrichment", "TSS enrichment")],
        summary_rows, empty_message="No TSS metrics available",
    )
    motif_rows = _top_motif_rows(data)
    expected = render_table(
        [("sample_id", "Sample"), ("expected_motif", "Expected motif"),
         ("expected_motif_status", "Expected motif status"), ("best_motif_id", "Best motif")],
        peak_rows, empty_message="No target motif metrics available",
    )
    expected += "<h3>Top AME motifs</h3>" + render_table(
        [("sample_id", "Sample"), ("rank", "Rank"), ("motif_id", "Motif"),
         ("motif_alt_id", "Alternate ID"), ("adjusted_p_value", "Adjusted p-value")],
        motif_rows, empty_message="No AME motifs available",
    )
    warnings = data.get("warnings", [])
    warning_rows = warnings if isinstance(warnings, list) else []
    warning_section = render_table(
        [("message", "Warning")], warning_rows, empty_message="No warnings recorded",
    )
    body = "".join((
        _section("run-overview", "Run overview", overview),
        '<p class="legend">Targets (solid); IgG controls (outlined) and hatched.</p>',
        _section("demultiplexing", "Demultiplexing", demux),
        _section("alignment", "Alignment and library QC", alignment),
        _section("peaks-frip", "Peaks and FRiP", peaks),
        _section("tss-enrichment", "TSS enrichment", tss),
        _section("motif-enrichment", "Motif enrichment", expected),
        _section("warnings", "Warnings", warning_section),
    ))
    return """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Consolidated QC dashboard</title><style>
body{font-family:system-ui,sans-serif;line-height:1.45;margin:0;color:#172033;background:#f8fafc}main{max-width:1100px;margin:auto;padding:1.5rem}section{background:#fff;border:1px solid #dbe3ee;border-radius:.5rem;padding:1rem;margin:1rem 0}h1,h2,h3{margin-top:0}.legend,.empty{color:#475569}table{border-collapse:collapse;width:100%;margin:.75rem 0}th,td{border:1px solid #dbe3ee;padding:.35rem;text-align:left;vertical-align:top}th{background:#eff6ff}.chart{width:100%;height:auto;background:#fff}.axis{stroke:#64748b}.chart-title{font-weight:700}
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
    for name in ("qc_dashboard.html", "qc_summary.tsv", "qc_summary.json", "top_motifs.tsv", "tss_profiles.tsv"):
        if not (outdir / name).is_file() or (outdir / name).stat().st_size == 0:
            raise DashboardInputError(f"dashboard output {name} is empty")


def write_outputs_atomically(data: Mapping[str, object], outdir: Path) -> None:
    """Publish a complete output set only after every artifact is non-empty."""
    outdir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{outdir.name}.tmp-", dir=outdir.parent))
    try:
        write_outputs(data, temporary)
        _validate_output_files(temporary)
        outdir.mkdir(parents=True, exist_ok=True)
        for name in ("qc_dashboard.html", "qc_summary.tsv", "qc_summary.json", "top_motifs.tsv", "tss_profiles.tsv"):
            os.replace(temporary / name, outdir / name)
    finally:
        if temporary.exists():
            for child in temporary.iterdir():
                child.unlink()
            temporary.rmdir()


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
    profiles: dict[str, list[tuple[int, float]]] = {}
    for path in _input_paths(directory, "*.tss_profile.tsv", label="TSS"):
        suffix = ".tss_profile.tsv"
        sample_id = path.name[:-len(suffix)]
        if not sample_id or sample_id in profiles:
            raise DashboardInputError(f"duplicate TSS profile sample_id {sample_id}")
        profiles[sample_id] = read_tss_profile(path)
    records = {sample_id: dict(record) for sample_id, record in statuses.items()}
    for sample_id, profile in profiles.items():
        record = records.setdefault(sample_id, {"sample_id": sample_id, "status": "computed"})
        record["profile"] = profile
        record["enrichment"] = calculate_tss_enrichment(profile)
    return dict(sorted(records.items()))


def _sample_id_for_ame(path: Path) -> str:
    """Recover the sample directory name from the AME producer's ``sample/ame`` tree."""
    return path.parent.parent.name if path.parent.name == "ame" else path.parent.name


def _read_top_motifs_from_directory(directory: Path) -> dict[str, list[dict[str, object]]]:
    records: dict[str, list[dict[str, object]]] = {}
    for path in _input_paths(directory, "ame.tsv", label="AME"):
        sample_id = _sample_id_for_ame(path)
        if not sample_id or sample_id in records:
            raise DashboardInputError(f"duplicate AME sample_id {sample_id}")
        records[sample_id] = read_top_ame(path)
    return dict(sorted(records.items()))


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
        _input_paths(args.insert_dir, "*.insert_size_distribution.tsv", label="insert")
        peaks = read_peak_metrics(_input_paths(args.peak_dir, "*.peak_qc.tsv", label="peak"))
        tss = _read_tss_metrics(args.tss_dir)
        motifs = _read_motif_metrics(
            _input_paths(args.motif_dir, "*.motif_qc.tsv", label="motif")
        )
        top_motifs = _read_top_motifs_from_directory(args.ame_dir)
        data = build_report_data(
            metadata, demultiplex, libraries, peaks, tss, motifs, top_motifs,
            annotation_status=args.annotation_status,
        )
        write_outputs_atomically(data, args.outdir)
    except DashboardInputError as error:
        print(f"qc_dashboard.py: error: {error}", file=__import__("sys").stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
