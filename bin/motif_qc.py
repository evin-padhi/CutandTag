#!/usr/bin/env python3
"""Summarize expected-motif AME enrichment and FIMO hit-position QC."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


SIGNIFICANCE_THRESHOLD = 0.05


class MotifQcError(ValueError):
    """Raised when MEME Suite tabular output violates the motif-QC contract."""


class AmeRecords(list["AmeRecord"]):
    """AME rows whose empty state means a parsed, thresholded AME report."""

    def __init__(self, records: Iterable["AmeRecord"] = (), *, complete_database: bool = False):
        super().__init__(records)
        self.complete_database = complete_database


@dataclass(frozen=True)
class AmeRecord:
    motif_id: str
    motif_alt_id: str
    adjusted_p_value: float
    p_value: float | None = None
    effect: float | None = None
    positive_sequences: int | None = None
    rank: int | None = None


@dataclass(frozen=True)
class FimoHit:
    motif_id: str
    motif_alt_id: str
    sequence_name: str
    start: int
    stop: int
    strand: str
    p_value: float | None = None
    q_value: float | None = None


def _normalise_header(value: str) -> str:
    # Preserve the distinction between AME's TP/FP counts and %TP/%FP rates.
    value = value.strip().lower().replace("%", "percent_")
    return re.sub(r"[^a-z0-9]+", "_", value).strip("_")


def _find_column(headers: list[str], aliases: Iterable[str]) -> str | None:
    wanted = set(aliases)
    return next((header for header in headers if header in wanted), None)


def _read_table(path: Path, *, label: str, required: tuple[set[str], ...]) -> list[dict[str, str]]:
    """Read a TSV table, accepting comments and MEME's commented header line."""

    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise MotifQcError(f"cannot read {label} file {path}: {error}") from error

    header: list[str] | None = None
    records: list[dict[str, str]] = []
    for line_number, raw_line in enumerate(lines, start=1):
        if not raw_line.strip():
            continue
        stripped = raw_line.lstrip()
        comment = stripped.startswith("#")
        candidate = stripped[1:].lstrip() if comment else raw_line
        fields = candidate.split("\t")
        normalised = [_normalise_header(field) for field in fields]
        is_header = all(_find_column(normalised, aliases) is not None for aliases in required)
        if header is None:
            if is_header:
                if len(set(normalised)) != len(normalised):
                    raise MotifQcError(f"{path}:{line_number}: duplicate {label} column names")
                header = normalised
            continue
        if comment:
            continue
        if is_header:
            # Some MEME versions include a prose commented header followed by a TSV header.
            header = normalised
            continue
        if len(fields) != len(header):
            raise MotifQcError(
                f"{path}:{line_number}: expected {len(header)} {label} columns, got {len(fields)}"
            )
        records.append(dict(zip(header, fields, strict=True)))
    if header is None:
        raise MotifQcError(f"{path}: no recognizable {label} TSV header")
    return records


def _ame_has_complete_database_marker(path: Path) -> bool:
    """Return pipeline-provided proof that AME was configured to report every motif.

    AME's standard TSV reports only enriched motifs.  The motif workflow adds this
    comment after it has configured an all-motifs report, allowing this utility to
    distinguish a regex absent from the database from a suppressed weak result.
    """

    try:
        return any(
            re.fullmatch(r"\s*#\s*motif_qc_complete_database\s*=\s*true\s*", line, re.IGNORECASE)
            for line in path.read_text(encoding="utf-8").splitlines()
        )
    except OSError as error:
        raise MotifQcError(f"cannot read AME file {path}: {error}") from error


def _text(record: dict[str, str], column: str | None, label: str) -> str:
    value = "" if column is None else record[column].strip()
    if not value:
        raise MotifQcError(f"{label} is required")
    return value


def _number(value: str, label: str, *, required: bool = False) -> float | None:
    value = value.strip()
    if not value or value == ".":
        if required:
            raise MotifQcError(f"{label} is required")
        return None
    try:
        parsed = float(value)
    except ValueError as error:
        raise MotifQcError(f"{label} must be numeric") from error
    if not math.isfinite(parsed):
        raise MotifQcError(f"{label} must be finite")
    return parsed


def _integer(value: str, label: str, *, required: bool = False) -> int | None:
    parsed = _number(value, label, required=required)
    if parsed is None:
        return None
    if not parsed.is_integer():
        raise MotifQcError(f"{label} must be an integer")
    return int(parsed)


def _probability(value: str, label: str, *, required: bool = False) -> float | None:
    parsed = _number(value, label, required=required)
    if parsed is not None and not 0 <= parsed <= 1:
        raise MotifQcError(f"{label} must be between 0 and 1")
    return parsed


def _nonnegative_number(value: str, label: str, *, required: bool = False) -> float | None:
    parsed = _number(value, label, required=required)
    if parsed is not None and parsed < 0:
        raise MotifQcError(f"{label} must be non-negative")
    return parsed


def read_ame(path: Path) -> list[AmeRecord]:
    """Read AME TSV output across common MEME Suite header variants."""

    motif_id_aliases = {"motif_id", "motif", "pattern_name"}
    # Modern AME reports a motif-level E-value rather than the removed
    # adj_p-value column.  It is AME's motif-level multiple-testing enrichment
    # statistic and is the comparable ranking/significance field here.
    adjusted_aliases = {"adj_p_value", "adjusted_p_value", "q_value", "e_value"}
    rows = _read_table(
        path,
        label="AME",
        required=(motif_id_aliases, adjusted_aliases),
    )
    complete_database = _ame_has_complete_database_marker(path)
    if not rows:
        # AME reports only significant motifs by default.  A header-only report
        # is therefore evidence of a completed AME run, not evidence of no peaks.
        return AmeRecords(complete_database=complete_database)
    headers = list(rows[0])
    motif_column = _find_column(headers, motif_id_aliases)
    adjusted_column = _find_column(headers, adjusted_aliases)
    alt_column = _find_column(headers, {"motif_alt_id", "motif_name", "alt_id"})
    p_column = _find_column(headers, {"p_value", "pvalue"})
    # AME's score is its primary effect statistic; enrichment is a useful fallback.
    effect_column = _find_column(headers, {"score", "enrichment", "effect", "odds_ratio"})
    positive_column = _find_column(
        headers, {"pos", "positive_sequences", "n_sequences", "num_sequences"}
    )
    negative_column = _find_column(headers, {"neg", "negative_sequences"})
    true_positive_column = _find_column(headers, {"tp", "true_positives"})
    false_positive_column = _find_column(headers, {"fp", "false_positives"})
    rank_column = _find_column(headers, {"rank"})
    parsed = AmeRecords(complete_database=complete_database)
    for row_number, row in enumerate(rows, start=1):
        positive = _integer(row[positive_column], f"AME row {row_number} positive sequences") if positive_column else None
        if positive is not None and positive < 0:
            raise MotifQcError(f"AME row {row_number} positive sequences must be non-negative")
        negative = _integer(
            row[negative_column],
            f"AME row {row_number} negative sequences",
        ) if negative_column else None
        if negative is not None and negative < 0:
            raise MotifQcError(
                f"AME row {row_number} negative sequences must be non-negative"
            )
        rank = _integer(row[rank_column], f"AME row {row_number} rank") if rank_column else None
        if rank is not None and rank < 1:
            raise MotifQcError(f"AME row {row_number} rank must be positive")
        if adjusted_column == "e_value":
            adjusted_value = _nonnegative_number(
                row[adjusted_column],
                f"AME row {row_number} E value",
                required=True,
            )
        else:
            adjusted_value = _probability(
                row[adjusted_column],
                f"AME row {row_number} adjusted P value",
                required=True,
            )
        assert adjusted_value is not None
        effect = (
            _number(row[effect_column], f"AME row {row_number} effect")
            if effect_column else None
        )
        if (
            effect is None
            and positive is not None
            and negative is not None
            and true_positive_column
            and false_positive_column
        ):
            true_positive = _integer(
                row[true_positive_column],
                f"AME row {row_number} true positives",
                required=True,
            )
            false_positive = _integer(
                row[false_positive_column],
                f"AME row {row_number} false positives",
                required=True,
            )
            assert true_positive is not None and false_positive is not None
            if (
                positive <= 0
                or negative <= 0
                or not 0 <= true_positive <= positive
                or not 0 <= false_positive <= negative
            ):
                raise MotifQcError(
                    f"AME row {row_number} Fisher counts are inconsistent "
                    "with positive/negative sequence totals"
                )
            positive_rate = true_positive / positive
            # A half-hit continuity correction keeps the ratio finite when
            # no control sequence is classified positive.
            negative_rate = (
                false_positive / negative
                if false_positive
                else 0.5 / negative
            )
            effect = positive_rate / negative_rate
        parsed.append(AmeRecord(
            motif_id=_text(row, motif_column, f"AME row {row_number} motif ID"),
            motif_alt_id=(row[alt_column].strip() if alt_column else ""),
            adjusted_p_value=adjusted_value,
            p_value=_probability(row[p_column], f"AME row {row_number} P value") if p_column else None,
            effect=effect,
            positive_sequences=positive,
            rank=rank,
        ))
    return parsed


def read_fimo(path: Path) -> list[FimoHit]:
    """Read FIMO TSV output across standard and human-readable header variants."""

    motif_id_aliases = {"motif_id", "pattern_name", "motif"}
    sequence_aliases = {"sequence_name", "sequence", "sequence_id"}
    rows = _read_table(
        path,
        label="FIMO",
        required=(motif_id_aliases, sequence_aliases, {"start"}, {"stop", "end"}),
    )
    if not rows:
        return []
    headers = list(rows[0])
    motif_column = _find_column(headers, motif_id_aliases)
    sequence_column = _find_column(headers, sequence_aliases)
    start_column = _find_column(headers, {"start"})
    stop_column = _find_column(headers, {"stop", "end"})
    alt_column = _find_column(headers, {"motif_alt_id", "pattern_accession", "motif_name"})
    strand_column = _find_column(headers, {"strand"})
    p_column = _find_column(headers, {"p_value", "pvalue"})
    q_column = _find_column(headers, {"q_value", "qvalue"})
    parsed: list[FimoHit] = []
    for row_number, row in enumerate(rows, start=1):
        start = _integer(row[start_column], f"FIMO row {row_number} start", required=True)
        stop = _integer(row[stop_column], f"FIMO row {row_number} stop", required=True)
        assert start is not None and stop is not None
        if start < 1 or stop < start:
            raise MotifQcError(f"FIMO row {row_number}: coordinates must be positive and inclusive")
        parsed.append(FimoHit(
            motif_id=_text(row, motif_column, f"FIMO row {row_number} motif ID"),
            motif_alt_id=(row[alt_column].strip() if alt_column else ""),
            sequence_name=_text(row, sequence_column, f"FIMO row {row_number} sequence name"),
            start=start,
            stop=stop,
            strand=(row[strand_column].strip() if strand_column else ""),
            p_value=_probability(row[p_column], f"FIMO row {row_number} P value") if p_column else None,
            q_value=_probability(row[q_column], f"FIMO row {row_number} Q value") if q_column else None,
        ))
    return parsed


def _matches(pattern: re.Pattern[str], record: AmeRecord) -> bool:
    return bool(pattern.search(record.motif_id) or pattern.search(record.motif_alt_id))


def summarize_expected_motif(
    ame_records: Iterable[AmeRecord],
    fimo_hits: Iterable[FimoHit],
    expected_motif: str,
    window: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Aggregate expected-motif enrichment, unique peak coverage, and hit positions."""

    if window <= 0:
        raise MotifQcError("--window must be a positive integer")
    try:
        pattern = re.compile(expected_motif)
    except re.error as error:
        raise MotifQcError(f"invalid expected-motif regular expression: {error}") from error
    parsed_ame_report = isinstance(ame_records, AmeRecords)
    complete_database = parsed_ame_report and ame_records.complete_database
    ame = list(ame_records)
    fimo = list(fimo_hits)
    peak_counts = {record.positive_sequences for record in ame if record.positive_sequences is not None}
    peak_count = next(iter(peak_counts)) if len(peak_counts) == 1 else None
    base = {
        "expected_motif": expected_motif,
        "ame_database_complete": complete_database,
        "peak_count": peak_count,
        "matching_motif_count": 0,
        "best_motif_id": None,
        "best_motif_alt_id": None,
        "best_motif_rank": None,
        "best_adjusted_p_value": None,
        "best_effect": None,
        "expected_hit_count": 0,
        "peaks_with_expected_motif": 0,
        "peak_fraction_with_expected_motif": None,
        "central_hit_count": 0,
        "central_hit_fraction": None,
    }
    if peak_count == 0:
        return ({"status": "no_peaks", **base}, [])
    if not ame:
        if parsed_ame_report:
            if complete_database:
                return ({"status": "motif_not_found", **base}, [])
            # AME's default report threshold suppresses non-significant motifs.
            return ({"status": "not_significant", **base}, [])
        return ({"status": "no_peaks", **base}, [])

    ranked = sorted(enumerate(ame), key=lambda item: (item[1].adjusted_p_value, item[1].motif_id, item[1].motif_alt_id, item[0]))
    matching = [(rank + 1, record) for rank, (_, record) in enumerate(ranked) if _matches(pattern, record)]
    if not matching:
        status = "motif_not_found" if complete_database else "not_significant"
        return ({"status": status, **base}, [])
    best_rank, best = matching[0]
    selected_ids = {record.motif_id for _, record in matching}
    selected_alts = {record.motif_alt_id for _, record in matching if record.motif_alt_id}
    selected_hits = [
        hit for hit in fimo
        if hit.motif_id in selected_ids or (hit.motif_alt_id and hit.motif_alt_id in selected_alts)
    ]
    center = (window + 1) / 2
    central_radius = window / 4
    positions: list[dict[str, object]] = []
    for hit in selected_hits:
        relative_position = (hit.start + hit.stop) / 2 - center
        positions.append({
            "motif_id": hit.motif_id,
            "motif_alt_id": hit.motif_alt_id,
            "sequence_name": hit.sequence_name,
            "start": hit.start,
            "stop": hit.stop,
            "strand": hit.strand,
            "p_value": hit.p_value,
            "q_value": hit.q_value,
            "relative_center_position": relative_position,
            "is_central": abs(relative_position) <= central_radius,
        })
    peaks_with_hits = len({hit.sequence_name for hit in selected_hits})
    central_hits = sum(position["is_central"] for position in positions)
    base.update({
        "matching_motif_count": len(matching),
        "best_motif_id": best.motif_id,
        "best_motif_alt_id": best.motif_alt_id or None,
        "best_motif_rank": best_rank,
        "best_adjusted_p_value": best.adjusted_p_value,
        "best_effect": best.effect,
        "expected_hit_count": len(selected_hits),
        "peaks_with_expected_motif": peaks_with_hits,
        "peak_fraction_with_expected_motif": (peaks_with_hits / peak_count) if peak_count else None,
        "central_hit_count": central_hits,
        "central_hit_fraction": central_hits / len(selected_hits) if selected_hits else 0.0,
    })
    status = "pass" if best.adjusted_p_value <= SIGNIFICANCE_THRESHOLD else "not_significant"
    return ({"status": status, **base}, positions)


def _output_paths(prefix: Path) -> dict[str, Path]:
    return {
        "json": Path(f"{prefix}.json"),
        "tsv": Path(f"{prefix}.tsv"),
        "positions": Path(f"{prefix}.motif_hit_positions.tsv"),
    }


def _format(value: object) -> str:
    if value is None:
        return ""
    if value is True:
        return "true"
    if value is False:
        return "false"
    return str(value)


def _backup_path(path: Path) -> Path:
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".backup", dir=path.parent)
    os.close(descriptor)
    backup = Path(name)
    backup.unlink()
    return backup


def _publish_outputs(staged: dict[Path, Path]) -> None:
    backups: dict[Path, Path] = {}
    published: list[Path] = []
    try:
        for final_path, staged_path in staged.items():
            if final_path.exists():
                backup = _backup_path(final_path)
                os.replace(final_path, backup)
                backups[final_path] = backup
            os.replace(staged_path, final_path)
            published.append(final_path)
    except BaseException:
        for final_path in published:
            final_path.unlink(missing_ok=True)
        for final_path, backup in backups.items():
            if backup.exists():
                os.replace(backup, final_path)
        raise
    else:
        for backup in backups.values():
            backup.unlink(missing_ok=True)


def write_motif_outputs(
    output_prefix: Path, summary: dict[str, object], positions: Iterable[dict[str, object]]
) -> None:
    """Write all expected-motif outputs as a rollback-safe publication set."""

    paths = _output_paths(output_prefix)
    position_rows = list(positions)
    metrics_tsv = ["metric\tvalue"]
    metrics_tsv.extend(f"{key}\t{_format(value)}" for key, value in summary.items())
    position_header = (
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tp_value\tq_value"
        "\trelative_center_position\tis_central"
    )
    position_tsv = [position_header]
    position_tsv.extend(
        "\t".join(_format(row[key]) for key in (
            "motif_id", "motif_alt_id", "sequence_name", "start", "stop", "strand",
            "p_value", "q_value", "relative_center_position", "is_central",
        ))
        for row in position_rows
    )
    rendered = {
        "json": json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        "tsv": "\n".join(metrics_tsv) + "\n",
        "positions": "\n".join(position_tsv) + "\n",
    }
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".motif-qc-", dir=output_prefix.parent) as stage:
        stage_path = Path(stage)
        staged: dict[Path, Path] = {}
        for key, final_path in paths.items():
            temporary_path = stage_path / final_path.name
            temporary_path.write_text(rendered[key], encoding="utf-8", newline="\n")
            staged[final_path] = temporary_path
        _publish_outputs(staged)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ame", required=True, type=Path, help="AME TSV result file")
    parser.add_argument("--fimo", required=True, type=Path, help="FIMO TSV result file")
    parser.add_argument("--expected-motif", required=True, help="Expected motif regular expression")
    parser.add_argument("--window", required=True, type=int, help="Peak-window width in bases")
    parser.add_argument("--output-prefix", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        summary, positions = summarize_expected_motif(
            read_ame(arguments.ame), read_fimo(arguments.fimo),
            arguments.expected_motif, arguments.window,
        )
        write_motif_outputs(arguments.output_prefix, summary, positions)
    except MotifQcError as error:
        parser.error(str(error))
    except OSError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
