#!/usr/bin/env python3
"""Parse and join the QC inputs used by the consolidated dashboard."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path
from typing import Mapping, Sequence

from motif_qc import read_ame


SCHEMA_VERSION = 1
TSS_BEFORE_BP = 3000
TSS_BIN_SIZE = 10
TSS_FLANK_BP = 100
TOP_MOTIF_LIMIT = 10


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
