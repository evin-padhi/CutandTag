#!/usr/bin/env python3
"""Summarize broad peaks and calculate fragment-based FRiP QC."""

from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


class PeakQcError(ValueError):
    """Raised when a BED/BEDPE input does not meet the QC contract."""


@dataclass(frozen=True)
class Peak:
    chrom: str
    start: int
    end: int
    name: str
    score: float | None = None
    signal_value: float | None = None

    @property
    def width(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class Fragment:
    chrom1: str
    start1: int
    end1: int
    chrom2: str
    start2: int
    end2: int
    name: str


def _validate_interval(chrom: str, start: int, end: int, label: str) -> None:
    if not chrom:
        raise PeakQcError(f"{label}: chromosome is required")
    if start < 0:
        raise PeakQcError(f"{label}: start must be non-negative")
    if end <= start:
        raise PeakQcError(f"{label}: interval must have positive width")


def _coerce_peak(record: Peak | Sequence[object], index: int) -> Peak:
    if isinstance(record, Peak):
        peak = Peak(
            record.chrom,
            record.start,
            record.end,
            record.name,
            _optional_number(record.score, f"peak {index}: score"),
            _optional_number(record.signal_value, f"peak {index}: signal value"),
        )
    else:
        if len(record) < 3:
            raise PeakQcError(f"peak {index}: expected at least three BED columns")
        try:
            peak = Peak(
                str(record[0]), int(record[1]), int(record[2]),
                str(record[3]) if len(record) > 3 else f"peak_{index}",
                _optional_number(record[4], f"peak {index}: score") if len(record) > 4 else None,
                _optional_number(record[6], f"peak {index}: signal value") if len(record) > 6 else None,
            )
        except PeakQcError:
            raise
        except (TypeError, ValueError) as error:
            raise PeakQcError(f"peak {index}: start and end must be integers") from error
    _validate_interval(peak.chrom, peak.start, peak.end, f"peak {index}")
    return peak


def _optional_number(value: object, label: str) -> float | None:
    if value is None or str(value) in {"", "."}:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise PeakQcError(f"{label} must be numeric") from error
    if not math.isfinite(parsed):
        raise PeakQcError(f"{label} must be finite")
    return parsed


def _coerce_fragment(record: Fragment | Sequence[object], index: int) -> Fragment:
    if isinstance(record, Fragment):
        fragment = record
    else:
        if len(record) < 7:
            raise PeakQcError(
                f"fragment {index}: BEDPE fragment name is required in column 7"
            )
        try:
            fragment = Fragment(
                str(record[0]), int(record[1]), int(record[2]),
                str(record[3]), int(record[4]), int(record[5]),
                _fragment_name(record[6], index),
            )
        except PeakQcError:
            raise
        except (TypeError, ValueError) as error:
            raise PeakQcError(f"fragment {index}: coordinates must be integers") from error
    _validate_interval(
        fragment.chrom1, fragment.start1, fragment.end1, f"fragment {index} mate 1"
    )
    _validate_interval(
        fragment.chrom2, fragment.start2, fragment.end2, f"fragment {index} mate 2"
    )
    name = _fragment_name(fragment.name, index)
    if name != fragment.name:
        return Fragment(
            fragment.chrom1,
            fragment.start1,
            fragment.end1,
            fragment.chrom2,
            fragment.start2,
            fragment.end2,
            name,
        )
    return fragment


def _fragment_name(value: object, index: int) -> str:
    name = "" if value is None else str(value).strip()
    if not name or name == ".":
        raise PeakQcError(f"fragment {index}: BEDPE fragment name is required in column 7")
    return name


def _normalise_peaks(peaks: Iterable[Peak | Sequence[object]]) -> list[Peak]:
    return [_coerce_peak(record, index) for index, record in enumerate(peaks, start=1)]


def _normalise_fragments(
    fragments: Iterable[Fragment | Sequence[object]],
) -> Iterator[Fragment]:
    return (
        _coerce_fragment(record, index) for index, record in enumerate(fragments, start=1)
    )


def _unique_fragments(fragments: Iterable[Fragment]) -> Iterator[Fragment]:
    """Keep one BEDPE record per required fragment name, preserving input order."""

    seen: set[str] = set()
    for fragment in fragments:
        if fragment.name not in seen:
            seen.add(fragment.name)
            yield fragment


def _quantile(values: list[int] | list[float], quantile: float) -> float | None:
    if not values:
        return None
    position = (len(values) - 1) * quantile
    lower = int(position)
    upper = min(lower + 1, len(values) - 1)
    fraction = position - lower
    return values[lower] + (values[upper] - values[lower]) * fraction


def _numeric_summary(values: list[float], prefix: str) -> dict[str, float | int | None]:
    sorted_values = sorted(values)
    return {
        f"{prefix}_count": len(sorted_values),
        f"{prefix}_min": sorted_values[0] if sorted_values else None,
        f"{prefix}_mean": sum(sorted_values) / len(sorted_values) if sorted_values else None,
        f"{prefix}_median": _quantile(sorted_values, 0.5),
        f"{prefix}_max": sorted_values[-1] if sorted_values else None,
        f"{prefix}_q25": _quantile(sorted_values, 0.25),
        f"{prefix}_q75": _quantile(sorted_values, 0.75),
    }


def _merged_peak_intervals(peaks: Iterable[Peak]) -> dict[str, list[tuple[int, int]]]:
    by_chromosome: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for peak in peaks:
        by_chromosome[peak.chrom].append((peak.start, peak.end))

    merged: dict[str, list[tuple[int, int]]] = {}
    for chrom, intervals in by_chromosome.items():
        combined: list[tuple[int, int]] = []
        for start, end in sorted(intervals):
            if combined and start <= combined[-1][1]:
                combined[-1] = (combined[-1][0], max(combined[-1][1], end))
            else:
                combined.append((start, end))
        merged[chrom] = combined
    return merged


def _fragment_intervals(fragment: Fragment) -> tuple[tuple[str, int, int], ...]:
    """Return the genomic span of a pair, preserving inter-chromosome mates."""

    if fragment.chrom1 == fragment.chrom2:
        return ((
            fragment.chrom1,
            min(fragment.start1, fragment.start2),
            max(fragment.end1, fragment.end2),
        ),)
    return (
        (fragment.chrom1, fragment.start1, fragment.end1),
        (fragment.chrom2, fragment.start2, fragment.end2),
    )


def _overlaps(start: int, end: int, other_start: int, other_end: int) -> bool:
    return start < other_end and end > other_start


@dataclass(frozen=True)
class _PeakIndex:
    starts: tuple[int, ...]
    prefix_max_ends: tuple[int, ...]
    intervals: tuple[tuple[int, int, int], ...]


def _index_peaks(peaks: list[Peak]) -> dict[str, _PeakIndex]:
    by_chromosome: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for index, peak in enumerate(peaks):
        by_chromosome[peak.chrom].append((peak.start, peak.end, index))

    indexed: dict[str, _PeakIndex] = {}
    for chrom, intervals in by_chromosome.items():
        ordered = tuple(sorted(intervals))
        maximum_end = -1
        prefix_max_ends: list[int] = []
        for _, end, _ in ordered:
            maximum_end = max(maximum_end, end)
            prefix_max_ends.append(maximum_end)
        indexed[chrom] = _PeakIndex(
            tuple(start for start, _, _ in ordered),
            tuple(prefix_max_ends),
            ordered,
        )
    return indexed


def _overlapping_peak_indexes(
    index: dict[str, _PeakIndex], chrom: str, start: int, end: int
) -> Iterable[int]:
    chrom_index = index.get(chrom)
    if chrom_index is None:
        return ()
    right = bisect_left(chrom_index.starts, end)
    left = bisect_right(chrom_index.prefix_max_ends, start, 0, right)
    return (
        peak_index
        for peak_start, peak_end, peak_index in chrom_index.intervals[left:right]
        if _overlaps(start, end, peak_start, peak_end)
    )


def _calculate_fragment_overlaps(
    fragments: Iterable[Fragment], peaks: list[Peak]
) -> tuple[dict[str, float | int], list[int]]:
    """Compute FRiP and every per-peak count in one fragment traversal."""

    index = _index_peaks(peaks)
    counts = [0] * len(peaks)
    total_fragments = 0
    fragments_in_peaks = 0
    for fragment in fragments:
        total_fragments += 1
        hit_peak_indexes: set[int] = set()
        for chrom, start, end in _fragment_intervals(fragment):
            hit_peak_indexes.update(_overlapping_peak_indexes(index, chrom, start, end))
        if hit_peak_indexes:
            fragments_in_peaks += 1
            for peak_index in hit_peak_indexes:
                counts[peak_index] += 1
    return (
        {
            "total_fragments": total_fragments,
            "fragments_in_peaks": fragments_in_peaks,
            "frip": (
                fragments_in_peaks / total_fragments if total_fragments else 0.0
            ),
        },
        counts,
    )


def summarize_peaks(peaks: Iterable[Peak | Sequence[object]]) -> dict[str, float | int | None]:
    """Return width statistics and union coverage for 0-based, half-open peaks."""

    normalised = _normalise_peaks(peaks)
    widths = sorted(peak.width for peak in normalised)
    covered_bases = sum(
        end - start
        for intervals in _merged_peak_intervals(normalised).values()
        for start, end in intervals
    )
    return {
        "peak_count": len(normalised),
        "total_covered_bases": covered_bases,
        "peak_width_min": widths[0] if widths else None,
        "peak_width_mean": sum(widths) / len(widths) if widths else None,
        "peak_width_median": _quantile(widths, 0.5),
        "peak_width_max": widths[-1] if widths else None,
        "peak_width_q25": _quantile(widths, 0.25),
        "peak_width_q75": _quantile(widths, 0.75),
        **_numeric_summary(
            [peak.score for peak in normalised if peak.score is not None],
            "peak_score",
        ),
        **_numeric_summary(
            [peak.signal_value for peak in normalised if peak.signal_value is not None],
            "signal_value",
        ),
    }


def calculate_frip(
    fragments: Iterable[Fragment | Sequence[object]],
    peaks: Iterable[Peak | Sequence[object]],
) -> dict[str, float | int]:
    """Calculate fragment-based FRiP, counting every BEDPE pair at most once."""

    normalised_fragments = _unique_fragments(_normalise_fragments(fragments))
    normalised_peaks = _normalise_peaks(peaks)
    return _calculate_fragment_overlaps(normalised_fragments, normalised_peaks)[0]


def read_peaks(path: Path) -> list[Peak]:
    """Read broadPeak/BED records, accepting comments and an empty file."""

    records: list[Peak] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise PeakQcError(f"cannot read peaks file {path}: {error}") from error
    for line_number, line in enumerate(lines, start=1):
        if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            continue
        fields = line.split("\t")
        try:
            records.append(_coerce_peak(fields, line_number))
        except PeakQcError as error:
            raise PeakQcError(f"{path}:{line_number}: {error}") from error
    return records


def read_fragments(path: Path) -> list[Fragment]:
    """Read name-collated BEDPE fragments, accepting comments and an empty file."""

    return list(iter_fragments(path))


def iter_fragments(path: Path) -> Iterator[Fragment]:
    """Yield validated BEDPE fragments without materializing the input file."""

    try:
        handle = path.open(encoding="utf-8")
    except OSError as error:
        raise PeakQcError(f"cannot read fragments file {path}: {error}") from error
    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if (
                not line
                or line.startswith("#")
                or line.startswith("track")
                or line.startswith("browser")
            ):
                continue
            fields = line.split("\t")
            try:
                yield _coerce_fragment(fields, line_number)
            except PeakQcError as error:
                raise PeakQcError(f"{path}:{line_number}: {error}") from error


def _output_paths(prefix: Path) -> dict[str, Path]:
    return {
        "json": Path(f"{prefix}.json"),
        "tsv": Path(f"{prefix}.tsv"),
        "histogram": Path(f"{prefix}.width_histogram.tsv"),
        "per_peak": Path(f"{prefix}.fragments_per_peak.tsv"),
    }


def _format_metric(value: float | int | None) -> str:
    return "" if value is None else str(value)


def _render_outputs(
    sample_id: str, peaks: list[Peak], fragments: Iterable[Fragment]
) -> dict[str, str]:
    summary = summarize_peaks(peaks)
    frip, counts = _calculate_fragment_overlaps(fragments, peaks)
    payload = {"sample_id": sample_id, **summary, **frip}
    metrics_tsv = ["metric\tvalue"]
    metrics_tsv.extend(f"{key}\t{_format_metric(value)}" for key, value in payload.items())
    histogram = Counter(peak.width for peak in peaks)
    histogram_tsv = ["width\tpeak_count"]
    histogram_tsv.extend(f"{width}\t{histogram[width]}" for width in sorted(histogram))
    per_peak_tsv = [
        "chrom\tstart\tend\tpeak_name\twidth\tscore\tsignal_value\tfragment_count"
    ]
    per_peak_tsv.extend(
        f"{peak.chrom}\t{peak.start}\t{peak.end}\t{peak.name}\t{peak.width}"
        f"\t{_format_metric(peak.score)}\t{_format_metric(peak.signal_value)}\t{count}"
        for peak, count in zip(peaks, counts, strict=True)
    )
    return {
        "json": json.dumps(payload, indent=2, sort_keys=True) + "\n",
        "tsv": "\n".join(metrics_tsv) + "\n",
        "histogram": "\n".join(histogram_tsv) + "\n",
        "per_peak": "\n".join(per_peak_tsv) + "\n",
    }


def _backup_path(path: Path) -> Path:
    descriptor, name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".backup", dir=path.parent
    )
    os.close(descriptor)
    backup = Path(name)
    backup.unlink()
    return backup


def _publish_outputs(staged: dict[Path, Path]) -> None:
    """Replace all output files, restoring prior files if publication fails."""

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


def write_qc_outputs(
    output_prefix: Path,
    sample_id: str,
    peaks: Iterable[Peak | Sequence[object]],
    fragments: Iterable[Fragment | Sequence[object]],
) -> None:
    """Write JSON and TSV QC outputs as one rollback-safe publication set."""

    paths = _output_paths(output_prefix)
    normalised_peaks = _normalise_peaks(peaks)
    normalised_fragments = _unique_fragments(_normalise_fragments(fragments))
    rendered = _render_outputs(sample_id, normalised_peaks, normalised_fragments)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".peak-qc-", dir=output_prefix.parent) as stage:
        stage_path = Path(stage)
        staged: dict[Path, Path] = {}
        for key, final_path in paths.items():
            temporary_path = stage_path / final_path.name
            temporary_path.write_text(rendered[key], encoding="utf-8", newline="\n")
            staged[final_path] = temporary_path
        _publish_outputs(staged)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fragments", required=True, type=Path, help="Name-collated BEDPE")
    parser.add_argument("--peaks", required=True, type=Path, help="Final broadPeak/BED file")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-prefix", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        peaks = read_peaks(arguments.peaks)
        fragments = iter_fragments(arguments.fragments)
        write_qc_outputs(arguments.output_prefix, arguments.sample_id, peaks, fragments)
    except PeakQcError as error:
        parser.error(str(error))
    except OSError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
