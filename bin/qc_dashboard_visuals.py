"""Deterministic visual styles and compact display formatting for QC dashboards."""

from dataclasses import dataclass
import hashlib
import math
from typing import Mapping, Sequence


TARGET_COLORS: Mapping[str, str] = {
    "IGG": "#8F8D87",
    "CTCF": "#2F78D1",
    "GATA1": "#1FAE7A",
    "RUNX1": "#F06432",
}

FALLBACK_COLORS = (
    "#7C3AED", "#0F766E", "#BE123C", "#A16207",
    "#0369A1", "#6D28D9", "#047857", "#B45309",
)

_TARGET_DASHES = ("solid", "4 3", "2 2", "8 3")
_CONTROL_DASH = "6 4"


@dataclass(frozen=True)
class SeriesStyle:
    color: str
    dash: str
    marker: str
    is_control: bool


def _normalized_target(assay_target: object) -> str:
    return str(assay_target).strip().upper()


def target_color(assay_target: object, *, is_control: bool = False) -> str:
    """Return a stable target color, reserving grey for controls."""
    if is_control:
        return TARGET_COLORS["IGG"]
    target = _normalized_target(assay_target)
    if target in TARGET_COLORS:
        return TARGET_COLORS[target]
    digest = hashlib.sha256(target.encode()).digest()
    return FALLBACK_COLORS[digest[0] % len(FALLBACK_COLORS)]


def series_style(sample_id: str, assay_target: object, *, is_control: bool) -> SeriesStyle:
    """Return a deterministic plotting style for one dashboard series."""
    if is_control:
        return SeriesStyle(target_color(assay_target, is_control=True), _CONTROL_DASH, "square", True)
    digest = hashlib.sha256(sample_id.encode()).digest()
    return SeriesStyle(
        target_color(assay_target),
        _TARGET_DASHES[digest[0] % len(_TARGET_DASHES)],
        "circle",
        False,
    )


def _three_significant(value: float) -> str:
    if value == 0:
        return "0"
    exponent = math.floor(math.log10(abs(value)))
    rounded = round(value, 2 - exponent)
    rounded_exponent = math.floor(math.log10(abs(rounded)))
    decimals = max(0, 2 - rounded_exponent)
    return f"{rounded:.{decimals}f}"


def _rounded_to_three_significant(value: float) -> float:
    if value == 0:
        return 0.0
    exponent = math.floor(math.log10(abs(value)))
    return round(value, 2 - exponent)


def format_significant(value: object, *, exact: bool = False, compact: bool = True) -> str:
    """Format a dashboard value with three significant digits where appropriate."""
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return str(value).lower()
    if exact and isinstance(value, (str, int)):
        return str(value)
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value) if exact else "NA"
    if not math.isfinite(numeric):
        return "NA"
    if exact:
        return str(value)
    absolute = abs(numeric)
    if absolute and absolute < 1e-4:
        return f"{numeric:.2e}"
    if compact:
        rounded = _rounded_to_three_significant(numeric)
        for threshold, suffix in ((1_000_000_000, "B"), (1_000_000, "M"), (1_000, "K")):
            if abs(rounded) >= threshold:
                return f"{_three_significant(rounded / threshold)}{suffix}"
    return _three_significant(numeric)


def _nonnegative_integer(value: object, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return value


def bin_weighted_series(
    series: Mapping[str, Sequence[tuple[int, int]]], *, bin_size: int = 250
) -> dict[str, list[dict[str, float | int]]]:
    """Aggregate every sample into shared bins from zero through the largest value."""
    if isinstance(bin_size, bool) or not isinstance(bin_size, int) or bin_size <= 0:
        raise ValueError("bin_size must be a positive integer")

    largest_value = None
    counts_by_sample: dict[str, dict[int, int]] = {}
    totals: dict[str, int] = {}
    for sample_id, observations in series.items():
        bin_counts: dict[int, int] = {}
        total = 0
        for position, count in observations:
            position = _nonnegative_integer(position, name="position")
            count = _nonnegative_integer(count, name="count")
            largest_value = position if largest_value is None else max(largest_value, position)
            bin_start = (position // bin_size) * bin_size
            bin_counts[bin_start] = bin_counts.get(bin_start, 0) + count
            total += count
        counts_by_sample[sample_id] = bin_counts
        totals[sample_id] = total

    if largest_value is None:
        return {sample_id: [] for sample_id in series}

    bin_starts = range(0, (largest_value // bin_size) * bin_size + bin_size, bin_size)
    return {
        sample_id: [
            {
                "bin_start": bin_start,
                "bin_end": bin_start + bin_size,
                "count": count,
                "percent": 100.0 * count / totals[sample_id] if totals[sample_id] else 0.0,
            }
            for bin_start in bin_starts
            for count in (counts_by_sample[sample_id].get(bin_start, 0),)
        ]
        for sample_id in series
    }


def histogram_ecdf(histogram: Sequence[tuple[int, int]]) -> list[dict[str, float | int]]:
    """Convert histogram counts to a sorted empirical cumulative distribution."""
    combined: dict[int, int] = {}
    for value, count in histogram:
        value = _nonnegative_integer(value, name="value")
        count = _nonnegative_integer(count, name="count")
        combined[value] = combined.get(value, 0) + count

    nonzero_items = [(value, count) for value, count in sorted(combined.items()) if count]
    total = sum(count for _, count in nonzero_items)
    cumulative = 0
    rows: list[dict[str, float | int]] = []
    for value, count in nonzero_items:
        cumulative += count
        rows.append({
            "value": value,
            "count": count,
            "cumulative_percent": 100.0 * cumulative / total,
        })
    return rows


def pack_endpoint_labels(
    endpoints: Sequence[tuple[str, float]], *, lower: float, upper: float, minimum_gap: float
) -> dict[str, float]:
    """Place endpoint labels with a minimum separation inside the plot bounds."""
    if not all(math.isfinite(value) for value in (lower, upper, minimum_gap)):
        raise ValueError("label bounds and minimum_gap must be finite")
    if lower > upper or minimum_gap < 0:
        raise ValueError("label bounds and minimum_gap are invalid")
    ordered = sorted(endpoints, key=lambda endpoint: endpoint[1])
    if minimum_gap * (len(ordered) - 1) > upper - lower:
        raise ValueError("impossible endpoint label layout")
    if not ordered:
        return {}

    positions = []
    for _, requested in ordered:
        if not math.isfinite(requested):
            raise ValueError("endpoint coordinates must be finite")
        position = min(max(requested, lower), upper)
        if positions:
            position = max(position, positions[-1] + minimum_gap)
        positions.append(position)

    overflow = positions[-1] - upper
    if overflow > 0:
        positions = [position - overflow for position in positions]
    for index in range(len(positions) - 2, -1, -1):
        positions[index] = min(positions[index], positions[index + 1] - minimum_gap)
    underflow = lower - positions[0]
    if underflow > 0:
        positions = [position + underflow for position in positions]
    return {label: position for (label, _), position in zip(ordered, positions)}
