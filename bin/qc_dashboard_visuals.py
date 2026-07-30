"""Deterministic visual styles and compact display formatting for QC dashboards."""

from dataclasses import dataclass
import hashlib
import html
import math
from typing import Callable, Mapping, Sequence


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


def _heatmap_fill(score: float) -> str:
    """Interpolate the heatmap's white-to-blue scale at a capped 0–60 score."""
    fraction = min(max(score, 0.0), 60.0) / 60.0
    start = (255, 255, 255)
    end = (29, 78, 216)
    channels = tuple(
        round(start_channel + fraction * (end_channel - start_channel))
        for start_channel, end_channel in zip(start, end)
    )
    return "#" + "".join(f"{channel:02X}" for channel in channels)


def render_motif_heatmap(matrix: Mapping[str, object]) -> str:
    """Render a target-only AME significance matrix as an accessible SVG."""
    raw_sample_ids = matrix.get("sample_ids", [])
    raw_motif_keys = matrix.get("motif_keys", [])
    raw_cells = matrix.get("cells", {})
    if (
        not isinstance(raw_sample_ids, Sequence)
        or isinstance(raw_sample_ids, (str, bytes))
        or not isinstance(raw_motif_keys, Sequence)
        or isinstance(raw_motif_keys, (str, bytes))
        or not isinstance(raw_cells, Mapping)
    ):
        raw_sample_ids, raw_motif_keys, raw_cells = [], [], {}

    sample_ids = [str(sample_id) for sample_id in raw_sample_ids]
    motif_keys = [
        (str(motif_key[0]), str(motif_key[1]))
        for motif_key in raw_motif_keys
        if (
            isinstance(motif_key, Sequence)
            and not isinstance(motif_key, (str, bytes))
            and len(motif_key) == 2
        )
    ]
    if not sample_ids or not motif_keys:
        return (
            '<div class="heatmap-scroll" role="region" '
            'aria-label="Motif enrichment heatmap" tabindex="0">'
            '<p class="empty" role="status">No motif enrichment data available.</p>'
            "</div>"
        )

    raw_targets = matrix.get("sample_targets", {})
    sample_targets = raw_targets if isinstance(raw_targets, Mapping) else {}
    cell_width = 64.0
    cell_height = 34.0
    row_label_width = 210.0
    top = 142.0
    right = 24.0
    legend_height = 66.0
    width = row_label_width + cell_width * len(sample_ids) + right
    height = top + cell_height * len(motif_keys) + legend_height

    sample_labels = []
    for column, sample_id in enumerate(sample_ids):
        x = row_label_width + (column + 0.5) * cell_width
        target = sample_targets.get(sample_id, sample_id.rsplit("_", 1)[-1])
        color = target_color(target)
        sample_labels.append(
            f'<text class="heatmap-sample-label" x="{x:.1f}" y="{top - 12:.1f}" '
            f'fill="{color}" text-anchor="start" '
            f'transform="rotate(-55 {x:.1f} {top - 12:.1f})">'
            f"{html.escape(sample_id)}</text>"
        )

    row_labels = []
    cells = []
    for row, motif_key in enumerate(motif_keys):
        motif_id, motif_alt_id = motif_key
        y = top + row * cell_height
        display_label = (
            f"{motif_alt_id} ({motif_id})" if motif_alt_id else motif_id
        )
        row_labels.append(
            f'<text class="heatmap-motif-label" x="{row_label_width - 10:.1f}" '
            f'y="{y + cell_height / 2 + 4:.1f}" text-anchor="end">'
            f"{html.escape(display_label)}</text>"
        )
        for column, sample_id in enumerate(sample_ids):
            raw_cell = raw_cells.get((sample_id, motif_key), {})
            cell = raw_cell if isinstance(raw_cell, Mapping) else {}
            numeric_score = _numeric(cell.get("score"))
            score = min(max(numeric_score or 0.0, 0.0), 60.0)
            adjusted = cell.get("adjusted_p_value")
            label = str(cell.get("label", "ns"))
            cognate = bool(cell.get("outlined"))
            x = row_label_width + column * cell_width
            adjusted_label = (
                _exact_number(float(adjusted))
                if _numeric(adjusted) is not None else "NA"
            )
            accessible_name = (
                f"{sample_id}; {motif_id} {motif_alt_id}".rstrip()
                + f"; adjusted p-value {adjusted_label}"
                + f"; −log10 adjusted p-value {_exact_number(score)}"
                + f"; {'cognate' if cognate else 'noncognate'}"
            )
            cell_class = "heatmap-cell cognate" if cognate else "heatmap-cell"
            text_color = "#FFFFFF" if score >= 36.0 else "#172033"
            cells.append(
                f'<g class="heatmap-cell-group" '
                f'aria-label="{html.escape(accessible_name, quote=True)}">'
                f"<title>{html.escape(accessible_name)}</title>"
                f'<rect class="{cell_class}" x="{x:.1f}" y="{y:.1f}" '
                f'width="{cell_width:.1f}" height="{cell_height:.1f}" '
                f'fill="{_heatmap_fill(score)}"/>'
                f'<text class="heatmap-cell-label" '
                f'x="{x + cell_width / 2:.1f}" '
                f'y="{y + cell_height / 2 + 4:.1f}" '
                f'fill="{text_color}" text-anchor="middle">'
                f"{html.escape(label, quote=False).replace('&gt;', '>')}</text>"
                "</g>"
            )

    legend_y = top + cell_height * len(motif_keys) + 22.0
    legend_x = row_label_width
    legend_width = min(240.0, cell_width * len(sample_ids))
    legend = (
        '<defs><linearGradient id="motif-significance-gradient">'
        '<stop offset="0%" stop-color="#FFFFFF"/>'
        '<stop offset="100%" stop-color="#1D4ED8"/>'
        "</linearGradient></defs>"
        f'<text class="heatmap-legend-title" x="{legend_x:.1f}" '
        f'y="{legend_y:.1f}">−log10 adjusted p-value (capped at 60)</text>'
        f'<rect class="heatmap-legend" x="{legend_x:.1f}" '
        f'y="{legend_y + 10:.1f}" width="{legend_width:.1f}" height="12" '
        'fill="url(#motif-significance-gradient)"/>'
        f'<text class="heatmap-legend-tick" x="{legend_x:.1f}" '
        f'y="{legend_y + 38:.1f}">0</text>'
        f'<text class="heatmap-legend-tick" '
        f'x="{legend_x + legend_width:.1f}" y="{legend_y + 38:.1f}" '
        'text-anchor="end">60</text>'
    )
    return (
        '<div class="heatmap-scroll" role="region" '
        'aria-label="Motif enrichment heatmap" tabindex="0">'
        f'<svg class="motif-heatmap" width="{width:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" role="img" '
        'aria-label="Motif enrichment heatmap" data-cell-width="64">'
        "<title>Motif enrichment heatmap</title>"
        + "".join(sample_labels)
        + "".join(row_labels)
        + "".join(cells)
        + legend
        + "</svg></div>"
    )


@dataclass(frozen=True)
class BarMetric:
    """Declarative configuration for one dashboard bar-chart metric."""

    key: str
    title: str
    axis_label: str
    value_multiplier: float = 1.0
    log10_axis: bool = False
    formatter: Callable[[object], str] = format_significant


def _numeric(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _exact_number(value: float) -> str:
    return format(value, ".15g")


def _axis_name(title: str, axis_label: str) -> str:
    return f"{title} ({axis_label})"


def _visible_number(value: float, axis_label: str) -> str:
    suffix = "%" if axis_label.strip() in {"%", "Percent"} else ""
    return f"{format_significant(value, compact=False)}{suffix}"


def _nice_linear_max(maximum: float) -> float:
    if maximum <= 0:
        return 1.0
    exponent = math.floor(math.log10(maximum))
    fraction = maximum / (10 ** exponent)
    nice_fraction = next(
        candidate for candidate in (1.0, 2.0, 5.0, 10.0)
        if candidate >= fraction
    )
    return nice_fraction * (10 ** exponent)


def _linear_ticks(maximum: float) -> list[float]:
    domain_max = _nice_linear_max(maximum)
    return [domain_max * index / 4.0 for index in range(5)]


def _svg_axis_y(
    ticks: Sequence[tuple[float, str]], *, left: float, right: float,
    plot_top: float, plot_bottom: float, domain_min: float, domain_max: float,
) -> str:
    span = domain_max - domain_min or 1.0
    fragments = []
    for value, label in ticks:
        y = plot_bottom - (value - domain_min) / span * (plot_bottom - plot_top)
        fragments.append(
            f'<line class="axis-grid" x1="{left:.1f}" y1="{y:.1f}" '
            f'x2="{right:.1f}" y2="{y:.1f}"/>'
            f'<text class="axis-tick-label axis-tick-y" x="{left - 8:.1f}" '
            f'y="{y + 4:.1f}" text-anchor="end">{html.escape(label)}</text>'
        )
    return "".join(fragments)


def _svg_axis_x_ticks(
    ticks: Sequence[tuple[float, str]], *, left: float, right: float,
    plot_bottom: float, domain_min: float, domain_max: float,
) -> str:
    span = domain_max - domain_min or 1.0
    fragments = []
    for value, label in ticks:
        x = left + (value - domain_min) / span * (right - left)
        fragments.append(
            f'<line class="axis-grid" x1="{x:.1f}" y1="32" '
            f'x2="{x:.1f}" y2="{plot_bottom:.1f}"/>'
            f'<text class="axis-tick-label axis-tick-x" x="{x:.1f}" '
            f'y="{plot_bottom + 18:.1f}" text-anchor="middle">{html.escape(label)}</text>'
        )
    return "".join(fragments)


def _empty_panel(title: str, axis_label: str) -> str:
    accessible_name = html.escape(_axis_name(title, axis_label), quote=True)
    return (
        f'<article class="qc-panel empty-panel" role="status" '
        f'aria-label="{accessible_name}"><h3>{html.escape(title)}</h3>'
        '<p class="empty">No numeric data available; values are NA.</p></article>'
    )


def _sample_label(
    sample_id: str, *, x: float, y: float, rotate: bool,
) -> str:
    escaped = html.escape(sample_id)
    if rotate:
        return (
            f'<text class="sample-label rotated" x="{x:.1f}" y="{y:.1f}" '
            f'transform="rotate(45 {x:.1f} {y:.1f})" '
            f'text-anchor="start">{escaped}</text>'
        )
    return (
        f'<text class="sample-label" x="{x:.1f}" y="{y:.1f}" '
        f'text-anchor="middle">{escaped}</text>'
    )


def _category_geometry(
    sample_ids: Sequence[str], *, minimum_slot_width: float
) -> tuple[bool, float, float, float]:
    """Size category slots and margins so rendered labels stay in the SVG."""
    longest = max((len(sample_id) for sample_id in sample_ids), default=0)
    estimated_text_width = longest * 7.0
    rotate_labels = longest > 14
    if not rotate_labels:
        return (
            False,
            max(minimum_slot_width, estimated_text_width + 18.0),
            18.0,
            42.0,
        )
    projected_text = estimated_text_width / math.sqrt(2.0)
    return (
        True,
        max(minimum_slot_width, projected_text + 18.0),
        projected_text + 18.0,
        projected_text + 34.0,
    )


def render_bar_panel(
    title: str,
    rows: Sequence[Mapping[str, object]],
    *,
    value_key: str,
    axis_label: str,
    value_multiplier: float = 1.0,
    log10_axis: bool = False,
) -> str:
    """Render one accessible, horizontally scalable QC bar panel."""
    if not rows:
        return _empty_panel(title, axis_label)
    prepared: list[tuple[Mapping[str, object], str, float | None]] = []
    for row in rows:
        number = _numeric(row.get(value_key))
        value = number * value_multiplier if number is not None else None
        if value is not None and (value < 0 or (log10_axis and value <= 0)):
            value = None
        prepared.append((row, str(row.get("sample_id", "")), value))
    numeric = [value for _, _, value in prepared if value is not None]
    if not numeric:
        return _empty_panel(title, axis_label)

    transformed = [math.log10(value) for value in numeric] if log10_axis else numeric
    if log10_axis:
        domain_min = math.floor(min(transformed)) - 1.0 if transformed else 0.0
        domain_max = math.ceil(max(transformed)) if transformed else 1.0
        tick_values = [
            domain_min + (domain_max - domain_min) * index / 4.0
            for index in range(5)
        ]
        ticks = [(value, f"10^{format(value, '.3g')}") for value in tick_values]
        displayed_axis_label = f"{axis_label} (log10)"
    else:
        domain_min = 0.0
        tick_values = _linear_ticks(max(numeric, default=0.0))
        domain_max = tick_values[-1]
        ticks = [
            (value, format_significant(value, compact=True))
            for value in tick_values
        ]
        displayed_axis_label = axis_label

    rotate_labels, slot_width, right_margin, label_space = _category_geometry(
        [sample_id for _, sample_id, _ in prepared],
        minimum_slot_width=56.0,
    )
    left = 62.0
    plot_top, plot_bottom = 32.0, 236.0
    height = plot_bottom + label_space
    width = max(500.0, left + right_margin + slot_width * len(prepared))
    plot_right = width - right_margin
    bar_width = slot_width * 0.62
    span = domain_max - domain_min or 1.0
    marks = []
    for index, (row, sample_id, value) in enumerate(prepared):
        center_x = left + slot_width * index + slot_width / 2.0
        label_y = plot_bottom + 18.0
        marks.append(
            _sample_label(
                sample_id, x=center_x, y=label_y, rotate=rotate_labels
            )
        )
        assay_target = str(row.get("assay_target", ""))
        is_control = bool(row.get("is_control"))
        style = series_style(sample_id, assay_target, is_control=is_control)
        escaped_sample = html.escape(sample_id, quote=True)
        escaped_target = html.escape(assay_target, quote=True)
        if value is None:
            marks.append(
                f'<text class="bar-value na-value" data-sample-id="{escaped_sample}" '
                f'x="{center_x:.1f}" y="{plot_bottom - 6:.1f}" '
                'text-anchor="middle">NA</text>'
            )
            continue
        plotted_value = math.log10(value) if log10_axis else value
        bar_height = max(
            0.0,
            (plotted_value - domain_min) / span * (plot_bottom - plot_top),
        )
        y = plot_bottom - bar_height
        exact = _exact_number(value)
        marks.append(
            f'<rect class="bar" data-sample-id="{escaped_sample}" '
            f'data-assay-target="{escaped_target}" '
            f'data-sample-kind="{"control" if is_control else "target"}" '
            f'x="{center_x - bar_width / 2.0:.1f}" y="{y:.1f}" '
            f'width="{bar_width:.1f}" height="{bar_height:.1f}" '
            f'fill="{style.color}" stroke="{style.color}" '
            f'stroke-width="{"3" if is_control else "1"}">'
            f'<title>{html.escape(sample_id)}: {exact}</title></rect>'
            f'<text class="bar-value" x="{center_x:.1f}" '
            f'y="{max(plot_top + 12.0, y - 5.0):.1f}" text-anchor="middle">'
            f'{html.escape(_visible_number(value, axis_label))}</text>'
        )

    accessible_name = _axis_name(title, axis_label)
    svg = (
        f'<svg class="panel-chart" width="{width:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" role="img" '
        f'aria-label="{html.escape(accessible_name, quote=True)}">'
        f'<title>{html.escape(accessible_name)}</title>'
        f'<text class="axis-title axis-title-y" x="16" '
        f'y="{(plot_top + plot_bottom) / 2:.1f}" '
        f'transform="rotate(-90 16 {(plot_top + plot_bottom) / 2:.1f})" '
        f'text-anchor="middle">{html.escape(displayed_axis_label)}</text>'
        + _svg_axis_y(
            ticks,
            left=left,
            right=plot_right,
            plot_top=plot_top,
            plot_bottom=plot_bottom,
            domain_min=domain_min,
            domain_max=domain_max,
        )
        + "".join(marks)
        + "</svg>"
    )
    return (
        f'<article class="qc-panel"><h3>{html.escape(title)}</h3>'
        f'<div class="panel-scroll" role="region" '
        f'aria-label="{html.escape(title, quote=True)} chart" tabindex="0">'
        f"{svg}</div></article>"
    )


def render_range_panel(
    title: str,
    rows: Sequence[Mapping[str, object]],
    *,
    minimum_key: str,
    q25_key: str,
    median_key: str,
    q75_key: str,
    maximum_key: str,
    axis_label: str,
) -> str:
    """Render minimum/maximum, interquartile range, and median per sample."""
    if not rows:
        return _empty_panel(title, axis_label)
    keys = (minimum_key, q25_key, median_key, q75_key, maximum_key)
    prepared = []
    for row in rows:
        values = tuple(_numeric(row.get(key)) for key in keys)
        valid = (
            all(value is not None and value >= 0 for value in values)
            and list(values) == sorted(values)
        )
        prepared.append((row, values if valid else None))
    observed = [
        value
        for _, values in prepared
        if values is not None
        for value in values
    ]
    if not observed:
        return _empty_panel(title, axis_label)
    domain_max = _nice_linear_max(max(observed, default=0.0))
    ticks = [
        (value, format_significant(value, compact=True))
        for value in _linear_ticks(domain_max)
    ]
    rotate_labels, slot_width, right_margin, label_space = _category_geometry(
        [str(row.get("sample_id", "")) for row, _ in prepared],
        minimum_slot_width=58.0,
    )
    left = 62.0
    plot_top, plot_bottom = 32.0, 236.0
    height = plot_bottom + label_space
    width = max(500.0, left + right_margin + slot_width * len(prepared))
    plot_right = width - right_margin
    marks = []
    for index, (row, values) in enumerate(prepared):
        sample_id = str(row.get("sample_id", ""))
        center_x = left + slot_width * index + slot_width / 2.0
        marks.append(
            _sample_label(
                sample_id, x=center_x, y=plot_bottom + 18.0,
                rotate=rotate_labels,
            )
        )
        if values is None:
            marks.append(
                f'<text class="na-value" x="{center_x:.1f}" '
                f'y="{plot_bottom - 6:.1f}" text-anchor="middle">NA</text>'
            )
            continue
        minimum, q25, median, q75, maximum = values

        def y(value: float) -> float:
            return plot_bottom - value / domain_max * (plot_bottom - plot_top)

        assay_target = str(row.get("assay_target", ""))
        style = series_style(
            sample_id, assay_target, is_control=bool(row.get("is_control"))
        )
        tooltip = (
            f"{sample_id}: min {_exact_number(minimum)}; "
            f"Q25 {_exact_number(q25)}; median {_exact_number(median)}; "
            f"Q75 {_exact_number(q75)}; max {_exact_number(maximum)}"
        )
        marks.append(
            f'<g data-assay-target="{html.escape(assay_target, quote=True)}">'
            f'<title>{html.escape(tooltip)}</title>'
            f'<line class="range-min-max" x1="{center_x:.1f}" y1="{y(minimum):.1f}" '
            f'x2="{center_x:.1f}" y2="{y(maximum):.1f}" stroke="{style.color}" '
            'stroke-width="2"/>'
            f'<line class="range-iqr" x1="{center_x:.1f}" y1="{y(q25):.1f}" '
            f'x2="{center_x:.1f}" y2="{y(q75):.1f}" stroke="{style.color}" '
            'stroke-width="8"/>'
            f'<circle class="range-median" cx="{center_x:.1f}" cy="{y(median):.1f}" '
            f'r="5" fill="{style.color}" stroke="#fff" stroke-width="1.5"/></g>'
        )
    accessible_name = _axis_name(title, axis_label)
    svg = (
        f'<svg class="panel-chart" width="{width:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" role="img" '
        f'aria-label="{html.escape(accessible_name, quote=True)}">'
        f'<title>{html.escape(accessible_name)}</title>'
        f'<text class="axis-title axis-title-y" x="16" '
        f'y="{(plot_top + plot_bottom) / 2:.1f}" '
        f'transform="rotate(-90 16 {(plot_top + plot_bottom) / 2:.1f})" '
        f'text-anchor="middle">{html.escape(axis_label)}</text>'
        + _svg_axis_y(
            ticks, left=left, right=plot_right, plot_top=plot_top,
            plot_bottom=plot_bottom, domain_min=0.0, domain_max=domain_max,
        )
        + "".join(marks)
        + "</svg>"
    )
    return (
        f'<article class="qc-panel"><h3>{html.escape(title)}</h3>'
        f'<div class="panel-scroll" role="region" '
        f'aria-label="{html.escape(title, quote=True)} chart" tabindex="0">'
        f"{svg}</div></article>"
    )


def render_scatter_panel(
    title: str,
    rows: Sequence[Mapping[str, object]],
    *,
    x_key: str,
    y_key: str,
    x_axis_label: str,
    y_axis_label: str,
    x_log10_axis: bool = False,
) -> str:
    """Render a labeled per-sample scatterplot with visible axes and leaders."""
    points = []
    for row in rows:
        x_value = _numeric(row.get(x_key))
        y_value = _numeric(row.get(y_key))
        if (
            x_value is None or y_value is None
            or x_value < 0 or y_value < 0
            or (x_log10_axis and x_value <= 0)
        ):
            continue
        points.append((row, math.log10(x_value) if x_log10_axis else x_value, y_value, x_value))
    displayed_x_axis = f"{x_axis_label} (log10)" if x_log10_axis else x_axis_label
    if not points:
        return _empty_panel(title, f"{displayed_x_axis}; {y_axis_label}")

    height = 330.0
    left, right, plot_top, plot_bottom = 66.0, 535.0, 32.0, 270.0
    label_x = right + 24.0
    longest_label_width = max(
        len(str(row.get("sample_id", ""))) * 7.0
        for row, _, _, _ in points
    )
    width = max(680.0, label_x + longest_label_width + 20.0)
    x_values = [point[1] for point in points]
    y_values = [point[2] for point in points]
    if x_log10_axis:
        x_min = math.floor(min(x_values))
        x_max = math.ceil(max(x_values))
        if x_min == x_max:
            x_min -= 1.0
        raw_x_ticks = [
            x_min + (x_max - x_min) * index / 4.0 for index in range(5)
        ]
        x_ticks = [(value, f"10^{format(value, '.3g')}") for value in raw_x_ticks]
    else:
        x_min = 0.0
        raw_x_ticks = _linear_ticks(max(x_values))
        x_max = raw_x_ticks[-1]
        x_ticks = [
            (value, format_significant(value, compact=True))
            for value in raw_x_ticks
        ]
    y_min = 0.0
    raw_y_ticks = _linear_ticks(max(y_values))
    y_max = raw_y_ticks[-1]
    y_ticks = [
        (value, format_significant(value, compact=True))
        for value in raw_y_ticks
    ]

    def x_position(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min or 1.0) * (right - left)

    def y_position(value: float) -> float:
        return plot_bottom - value / (y_max or 1.0) * (plot_bottom - plot_top)

    requested_labels = [
        (str(row.get("sample_id", "")), y_position(y_value))
        for row, _, y_value, _ in points
    ]
    minimum_gap = min(16.0, (plot_bottom - plot_top) / max(1, len(points) - 1))
    label_positions = pack_endpoint_labels(
        requested_labels,
        lower=plot_top + 4.0,
        upper=plot_bottom - 4.0,
        minimum_gap=minimum_gap,
    )
    marks = []
    for row, plotted_x, y_value, original_x in points:
        sample_id = str(row.get("sample_id", ""))
        assay_target = str(row.get("assay_target", ""))
        style = series_style(
            sample_id, assay_target, is_control=bool(row.get("is_control"))
        )
        point_x, point_y = x_position(plotted_x), y_position(y_value)
        label_y = label_positions[sample_id]
        tooltip = (
            f"{sample_id}: x {_exact_number(original_x)}; "
            f"y {_exact_number(y_value)}"
        )
        marks.append(
            f'<line class="scatter-leader" x1="{point_x:.1f}" y1="{point_y:.1f}" '
            f'x2="{label_x - 4.0:.1f}" y2="{label_y:.1f}" '
            f'stroke="{style.color}" stroke-width="1"/>'
            f'<circle class="scatter-point" data-assay-target="'
            f'{html.escape(assay_target, quote=True)}" cx="{point_x:.1f}" '
            f'cy="{point_y:.1f}" r="5" fill="{style.color}">'
            f'<title>{html.escape(tooltip)}</title></circle>'
            f'<text class="scatter-label" x="{label_x:.1f}" y="{label_y + 4.0:.1f}" '
            f'fill="{style.color}">{html.escape(sample_id)}</text>'
        )
    accessible_name = f"{title}: {displayed_x_axis} by {y_axis_label}"
    svg = (
        f'<svg class="panel-chart scatter-chart" width="{width:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" '
        f'role="img" aria-label="{html.escape(accessible_name, quote=True)}">'
        f'<title>{html.escape(accessible_name)}</title>'
        f'<text class="axis-title axis-title-x" x="{(left + right) / 2:.1f}" '
        f'y="{height - 8:.1f}" text-anchor="middle">'
        f'{html.escape(displayed_x_axis)}</text>'
        f'<text class="axis-title axis-title-y" x="16" '
        f'y="{(plot_top + plot_bottom) / 2:.1f}" '
        f'transform="rotate(-90 16 {(plot_top + plot_bottom) / 2:.1f})" '
        f'text-anchor="middle">{html.escape(y_axis_label)}</text>'
        + _svg_axis_x_ticks(
            x_ticks, left=left, right=right, plot_bottom=plot_bottom,
            domain_min=x_min, domain_max=x_max,
        )
        + _svg_axis_y(
            y_ticks, left=left, right=right, plot_top=plot_top,
            plot_bottom=plot_bottom, domain_min=y_min, domain_max=y_max,
        )
        + "".join(marks)
        + "</svg>"
    )
    return (
        f'<article class="qc-panel"><h3>{html.escape(title)}</h3>'
        f'<div class="panel-scroll" role="region" '
        f'aria-label="{html.escape(title, quote=True)} chart" tabindex="0">'
        f"{svg}</div></article>"
    )


def _series_metadata(
    sample_id: str, metadata: Mapping[str, Mapping[str, object]]
) -> tuple[str, bool, SeriesStyle]:
    sample_metadata = metadata.get(sample_id, {})
    assay_target = str(sample_metadata.get("assay_target", ""))
    is_control = bool(sample_metadata.get("is_control"))
    return (
        assay_target,
        is_control,
        series_style(sample_id, assay_target, is_control=is_control),
    )


def _allocate_series_styles(
    sample_ids: Sequence[str],
    metadata: Mapping[str, Mapping[str, object]],
) -> dict[str, SeriesStyle]:
    """Allocate collision-free line styles within each assay target."""
    grouped: dict[str, list[str]] = {}
    for sample_id in sample_ids:
        assay_target, _, _ = _series_metadata(sample_id, metadata)
        grouped.setdefault(_normalized_target(assay_target), []).append(sample_id)

    allocated: dict[str, SeriesStyle] = {}
    for target_sample_ids in grouped.values():
        for index, sample_id in enumerate(sorted(target_sample_ids)):
            _, _, base_style = _series_metadata(sample_id, metadata)
            dash = "none" if index == 0 else f"{index + 1} {index + 2}"
            allocated[sample_id] = SeriesStyle(
                color=base_style.color,
                dash=dash,
                marker=base_style.marker,
                is_control=base_style.is_control,
            )
    return allocated


def _line_marker(
    style: SeriesStyle, *, x: float, y: float, tooltip: str,
    css_class: str = "series-point", attributes: str = "",
) -> str:
    escaped_tooltip = html.escape(tooltip)
    if style.marker == "square":
        return (
            f'<rect class="{css_class}" {attributes}x="{x - 3.0:.1f}" '
            f'y="{y - 3.0:.1f}" width="6" height="6" fill="{style.color}">'
            f"<title>{escaped_tooltip}</title></rect>"
        )
    return (
        f'<circle class="{css_class}" {attributes}cx="{x:.1f}" cy="{y:.1f}" '
        f'r="3" fill="{style.color}"><title>{escaped_tooltip}</title></circle>'
    )


def _direct_label_geometry(
    endpoints: Sequence[tuple[str, float]], *, plot_top: float,
    minimum_plot_bottom: float, minimum_gap: float = 14.0,
) -> tuple[float, dict[str, float]]:
    """Pack endpoint labels, expanding the plot rather than omitting a sample."""
    required_span = minimum_gap * max(0, len(endpoints) - 1)
    plot_bottom = max(minimum_plot_bottom, plot_top + required_span + 8.0)
    positions = pack_endpoint_labels(
        endpoints,
        lower=plot_top + 4.0,
        upper=plot_bottom - 4.0,
        minimum_gap=minimum_gap,
    )
    return plot_bottom, positions


def _endpoint_label(
    sample_id: str, *, endpoint_x: float, endpoint_y: float,
    label_x: float, label_y: float, style: SeriesStyle,
    attributes: str = "",
) -> str:
    leader = ""
    if abs(endpoint_y - label_y) >= 0.5:
        leader = (
            f'<line class="endpoint-leader" x1="{endpoint_x:.1f}" '
            f'y1="{endpoint_y:.1f}" x2="{label_x - 4.0:.1f}" '
            f'y2="{label_y:.1f}" stroke="{style.color}" stroke-width="1"/>'
        )
    return (
        leader
        + f'<text class="endpoint-label" data-sample-id="'
        f'{html.escape(sample_id, quote=True)}" {attributes}'
        f'x="{label_x:.1f}" '
        f'y="{label_y:.1f}" fill="{style.color}" dominant-baseline="middle">'
        f"{html.escape(sample_id)}</text>"
    )


def _line_legend_entry(
    sample_id: str, assay_target: str, style: SeriesStyle
) -> str:
    kind = "control" if style.is_control else "target"
    return (
        f'<li data-sample-kind="{kind}">'
        '<svg class="series-swatch" viewBox="0 0 36 10" '
        'aria-hidden="true" focusable="false">'
        f'<line x1="0" y1="5" x2="36" y2="5" stroke="{style.color}" '
        f'stroke-width="2" stroke-dasharray="{style.dash}"/></svg>'
        f'<span>{html.escape(sample_id)} ({html.escape(assay_target)})</span></li>'
    )


def render_binned_distribution(
    title: str,
    series: Mapping[str, Sequence[Mapping[str, object]]],
    metadata: Mapping[str, Mapping[str, object]],
    *,
    x_axis_label: str,
) -> str:
    """Render cohort-normalized binned distributions with direct sample labels."""
    prepared: dict[str, list[tuple[float, float, int, int, int]]] = {}
    for sample_id, rows in series.items():
        points = []
        for row in rows:
            bin_start = _numeric(row.get("bin_start"))
            bin_end = _numeric(row.get("bin_end"))
            percent = _numeric(row.get("percent"))
            count = _numeric(row.get("count"))
            if (
                bin_start is None or bin_end is None or percent is None
                or count is None or bin_start < 0 or bin_end <= bin_start
                or percent < 0 or count < 0
            ):
                continue
            points.append((
                (bin_start + bin_end) / 2.0,
                percent,
                int(bin_start),
                int(bin_end),
                int(count),
            ))
        if points:
            prepared[str(sample_id)] = sorted(points)
    if not prepared:
        return _empty_panel(title, f"{x_axis_label}; Percent")

    all_points = [point for points in prepared.values() for point in points]
    x_max = max(point[3] for point in all_points)
    y_max = _nice_linear_max(max(point[1] for point in all_points))
    left, plot_top, minimum_plot_bottom = 62.0, 32.0, 236.0
    right, label_x = 535.0, 557.0

    def provisional_y(value: float) -> float:
        return minimum_plot_bottom - value / y_max * (
            minimum_plot_bottom - plot_top
        )

    anchor_indexes = {
        sample_id: next(
            (
                index
                for index in range(len(points) - 1, -1, -1)
                if points[index][4] > 0
            ),
            len(points) - 1,
        )
        for sample_id, points in prepared.items()
    }
    endpoint_requests = [
        (sample_id, provisional_y(points[anchor_indexes[sample_id]][1]))
        for sample_id, points in prepared.items()
    ]
    plot_bottom, label_positions = _direct_label_geometry(
        endpoint_requests,
        plot_top=plot_top,
        minimum_plot_bottom=minimum_plot_bottom,
    )
    width = max(
        680.0,
        label_x + max(len(sample_id) for sample_id in prepared) * 7.0 + 18.0,
    )
    height = plot_bottom + 48.0

    def x_position(value: float) -> float:
        return left + value / (x_max or 1.0) * (right - left)

    def y_position(value: float) -> float:
        return plot_bottom - value / y_max * (plot_bottom - plot_top)

    x_ticks_raw = [x_max * index / 4.0 for index in range(5)]
    x_ticks = [
        (value, format_significant(value, compact=True))
        for value in x_ticks_raw
    ]
    percent_axis_label = (
        "Percent of read pairs"
        if "insert" in title.lower()
        else "Percent of peaks"
    )
    y_ticks = [
        (value, format_significant(value, compact=False))
        for value in _linear_ticks(y_max)
    ]
    marks = []
    legend = []
    styles = _allocate_series_styles(list(prepared), metadata)
    for sample_id, points in sorted(prepared.items()):
        assay_target, _, _ = _series_metadata(sample_id, metadata)
        style = styles[sample_id]
        plotted = [
            (x_position(midpoint), y_position(percent))
            for midpoint, percent, _, _, _ in points
        ]
        attributes = (
            f'data-sample-id="{html.escape(sample_id, quote=True)}" '
            f'data-assay-target="{html.escape(assay_target, quote=True)}" '
            f'data-sample-kind="{"control" if style.is_control else "target"}" '
            f'data-marker="{style.marker}" '
        )
        marks.append(
            f'<polyline class="distribution-trace" {attributes}'
            f'points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in plotted)}" '
            f'fill="none" stroke="{style.color}" stroke-width="2" '
            f'stroke-dasharray="{style.dash}"/>'
        )
        for (midpoint, percent, bin_start, bin_end, count), (x, y) in zip(
            points, plotted
        ):
            tooltip = (
                f"{sample_id}: [{bin_start}, {bin_end}) bp; "
                f"{_exact_number(percent)}%; count {count}"
            )
            point_attributes = (
                f'data-sample-id="{html.escape(sample_id, quote=True)}" '
                f'data-percent="{_exact_number(percent)}" '
                f'data-bin-midpoint="{_exact_number(midpoint)}" '
            )
            marks.append(
                _line_marker(
                    style, x=x, y=y, tooltip=tooltip,
                    attributes=point_attributes,
                )
            )
        anchor_index = anchor_indexes[sample_id]
        endpoint_x, endpoint_y = plotted[anchor_index]
        anchor_midpoint = points[anchor_index][0]
        marks.append(
            _endpoint_label(
                sample_id,
                endpoint_x=endpoint_x,
                endpoint_y=endpoint_y,
                label_x=label_x,
                label_y=label_positions[sample_id],
                style=style,
                attributes=(
                    f'data-anchor-bin-midpoint="{_exact_number(anchor_midpoint)}" '
                ),
            )
        )
        legend.append(_line_legend_entry(sample_id, assay_target, style))

    accessible_name = f"{title}: {x_axis_label} by percent"
    svg = (
        f'<svg class="panel-chart distribution-chart" width="{width:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" role="img" '
        f'aria-label="{html.escape(accessible_name, quote=True)}">'
        f"<title>{html.escape(accessible_name)}</title>"
        f'<text class="axis-title axis-title-x" x="{(left + right) / 2:.1f}" '
        f'y="{height - 8:.1f}" text-anchor="middle">'
        f"{html.escape(x_axis_label)}</text>"
        f'<text class="axis-title axis-title-y" x="16" '
        f'y="{(plot_top + plot_bottom) / 2:.1f}" '
        f'transform="rotate(-90 16 {(plot_top + plot_bottom) / 2:.1f})" '
        f'text-anchor="middle">{percent_axis_label}</text>'
        + _svg_axis_x_ticks(
            x_ticks, left=left, right=right, plot_bottom=plot_bottom,
            domain_min=0.0, domain_max=x_max,
        )
        + _svg_axis_y(
            y_ticks, left=left, right=right, plot_top=plot_top,
            plot_bottom=plot_bottom, domain_min=0.0, domain_max=y_max,
        )
        + "".join(marks)
        + "</svg>"
    )
    return (
        f'<article class="qc-panel"><h3>{html.escape(title)}</h3>'
        f'<div class="panel-scroll" role="region" '
        f'aria-label="{html.escape(title, quote=True)} chart" tabindex="0">'
        f"{svg}</div><ol class=\"series-legend compact\">"
        + "".join(legend)
        + "</ol></article>"
    )


def render_ecdf(
    title: str,
    series: Mapping[str, Sequence[Mapping[str, object]]],
    metadata: Mapping[str, Mapping[str, object]],
    *,
    x_axis_label: str,
    zero_origin: bool,
) -> str:
    """Render ECDF traces, reserving a separate x-axis slot for exact zero."""
    prepared: dict[str, list[tuple[int, float, int]]] = {}
    for sample_id, rows in series.items():
        points = []
        for row in rows:
            value = _numeric(row.get("value"))
            cumulative = _numeric(row.get("cumulative_percent"))
            count = _numeric(row.get("count"))
            if (
                value is None or cumulative is None or count is None
                or value < 0 or cumulative < 0 or count < 0
                or int(value) != value
            ):
                continue
            points.append((int(value), cumulative, int(count)))
        if points:
            prepared[str(sample_id)] = sorted(points)
    if not prepared:
        return _empty_panel(title, f"{x_axis_label}; Cumulative percent")

    maximum_value = max(
        value for points in prepared.values() for value, _, _ in points
    )
    maximum_power = max(0, math.ceil(math.log10(maximum_value))) if maximum_value else 0
    left, plot_top, minimum_plot_bottom = 62.0, 32.0, 236.0
    right, label_x = 535.0, 557.0
    x_domain_max = float(maximum_power + 1 if zero_origin else maximum_power)
    if x_domain_max <= 0:
        x_domain_max = 1.0

    def transformed_x(value: int) -> float:
        if zero_origin:
            return 0.0 if value == 0 else 1.0 + math.log10(value)
        return math.log10(max(value, 1))

    endpoint_requests = [
        (
            sample_id,
            minimum_plot_bottom - points[-1][1] / 100.0
            * (minimum_plot_bottom - plot_top),
        )
        for sample_id, points in prepared.items()
    ]
    plot_bottom, label_positions = _direct_label_geometry(
        endpoint_requests,
        plot_top=plot_top,
        minimum_plot_bottom=minimum_plot_bottom,
    )
    width = max(
        680.0,
        label_x + max(len(sample_id) for sample_id in prepared) * 7.0 + 18.0,
    )
    height = plot_bottom + 48.0

    def x_position(value: int) -> float:
        return left + transformed_x(value) / x_domain_max * (right - left)

    def y_position(value: float) -> float:
        return plot_bottom - value / 100.0 * (plot_bottom - plot_top)

    x_ticks = []
    if zero_origin:
        x_ticks.append((0.0, "0"))
        x_ticks.extend(
            (float(power + 1), str(10 ** power))
            for power in range(maximum_power + 1)
        )
    else:
        x_ticks.extend(
            (float(power), str(10 ** power))
            for power in range(maximum_power + 1)
        )
    y_ticks = [(value, format(value, ".0f")) for value in (0, 25, 50, 75, 100)]
    marks = []
    legend = []
    styles = _allocate_series_styles(list(prepared), metadata)
    for sample_id, points in sorted(prepared.items()):
        assay_target, _, _ = _series_metadata(sample_id, metadata)
        style = styles[sample_id]
        plotted = [
            (x_position(value), y_position(cumulative))
            for value, cumulative, _ in points
        ]
        attributes = (
            f'data-sample-id="{html.escape(sample_id, quote=True)}" '
            f'data-assay-target="{html.escape(assay_target, quote=True)}" '
            f'data-sample-kind="{"control" if style.is_control else "target"}" '
            f'data-marker="{style.marker}" '
        )
        marks.append(
            f'<polyline class="ecdf-trace" {attributes}'
            f'points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in plotted)}" '
            f'fill="none" stroke="{style.color}" stroke-width="2" '
            f'stroke-dasharray="{style.dash}"/>'
        )
        for (value, cumulative, count), (x, y) in zip(points, plotted):
            tooltip = (
                f"{sample_id}: {value} fragments; "
                f"{_exact_number(cumulative)}% cumulative; count {count}"
            )
            marks.append(
                _line_marker(style, x=x, y=y, tooltip=tooltip)
            )
        endpoint_x, endpoint_y = plotted[-1]
        marks.append(
            _endpoint_label(
                sample_id,
                endpoint_x=endpoint_x,
                endpoint_y=endpoint_y,
                label_x=label_x,
                label_y=label_positions[sample_id],
                style=style,
            )
        )
        legend.append(_line_legend_entry(sample_id, assay_target, style))

    accessible_name = f"{title}: {x_axis_label} ECDF"
    svg = (
        f'<svg class="panel-chart ecdf-chart" data-zero-origin="'
        f'{"true" if zero_origin else "false"}" width="{width:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" role="img" '
        f'aria-label="{html.escape(accessible_name, quote=True)}">'
        f"<title>{html.escape(accessible_name)}</title>"
        f'<text class="axis-title axis-title-x" x="{(left + right) / 2:.1f}" '
        f'y="{height - 8:.1f}" text-anchor="middle">'
        f"{html.escape(x_axis_label)} (log10 positive values)</text>"
        f'<text class="axis-title axis-title-y" x="16" '
        f'y="{(plot_top + plot_bottom) / 2:.1f}" '
        f'transform="rotate(-90 16 {(plot_top + plot_bottom) / 2:.1f})" '
        'text-anchor="middle">Cumulative percent of peaks</text>'
        + _svg_axis_x_ticks(
            x_ticks, left=left, right=right, plot_bottom=plot_bottom,
            domain_min=0.0, domain_max=x_domain_max,
        )
        + _svg_axis_y(
            y_ticks, left=left, right=right, plot_top=plot_top,
            plot_bottom=plot_bottom, domain_min=0.0, domain_max=100.0,
        )
        + "".join(marks)
        + "</svg>"
    )
    return (
        f'<article class="qc-panel"><h3>{html.escape(title)}</h3>'
        f'<div class="panel-scroll" role="region" '
        f'aria-label="{html.escape(title, quote=True)} chart" tabindex="0">'
        f"{svg}</div><ol class=\"series-legend compact\">"
        + "".join(legend)
        + "</ol></article>"
    )


def render_profile_chart(
    title: str,
    series: Mapping[str, Sequence[tuple[int, float | None]]],
    metadata: Mapping[str, Mapping[str, object]],
    *,
    x_axis_label: str,
    y_axis_label: str,
) -> str:
    """Render target-colored profiles with a TSS reference and direct labels."""
    prepared: dict[str, list[tuple[float, float]]] = {}
    for sample_id, rows in series.items():
        points = []
        for position, signal in rows:
            numeric_position = _numeric(position)
            numeric_signal = _numeric(signal)
            if (
                numeric_position is None or numeric_signal is None
                or numeric_signal < 0
            ):
                continue
            points.append((numeric_position, numeric_signal))
        if points:
            prepared[str(sample_id)] = sorted(points)
    if not prepared:
        return _empty_panel(title, f"{x_axis_label}; {y_axis_label}")

    all_points = [point for points in prepared.values() for point in points]
    x_min = min(point[0] for point in all_points)
    x_max = max(point[0] for point in all_points)
    y_max = _nice_linear_max(max(point[1] for point in all_points))
    left, plot_top, minimum_plot_bottom = 62.0, 32.0, 236.0
    right, label_x = 535.0, 557.0

    def provisional_y(value: float) -> float:
        return minimum_plot_bottom - value / y_max * (
            minimum_plot_bottom - plot_top
        )

    plot_bottom, label_positions = _direct_label_geometry(
        [
            (sample_id, provisional_y(points[-1][1]))
            for sample_id, points in prepared.items()
        ],
        plot_top=plot_top,
        minimum_plot_bottom=minimum_plot_bottom,
    )
    width = max(
        680.0,
        label_x + max(len(sample_id) for sample_id in prepared) * 7.0 + 18.0,
    )
    height = plot_bottom + 48.0

    def x_position(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min or 1.0) * (right - left)

    def y_position(value: float) -> float:
        return plot_bottom - value / y_max * (plot_bottom - plot_top)

    x_ticks_raw = [
        x_min + (x_max - x_min) * index / 4.0 for index in range(5)
    ]
    x_ticks = [
        (value, format_significant(value, compact=True))
        for value in x_ticks_raw
    ]
    y_ticks = [
        (value, format_significant(value, compact=True))
        for value in _linear_ticks(y_max)
    ]
    marks = []
    legend = []
    styles = _allocate_series_styles(list(prepared), metadata)
    for sample_id, points in sorted(prepared.items()):
        assay_target, _, _ = _series_metadata(sample_id, metadata)
        style = styles[sample_id]
        plotted = [
            (x_position(position), y_position(signal))
            for position, signal in points
        ]
        attributes = (
            f'data-sample-id="{html.escape(sample_id, quote=True)}" '
            f'data-assay-target="{html.escape(assay_target, quote=True)}" '
            f'data-sample-kind="{"control" if style.is_control else "target"}" '
            f'data-marker="{style.marker}" '
        )
        marks.append(
            f'<polyline class="profile-trace" {attributes}'
            f'points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in plotted)}" '
            f'fill="none" stroke="{style.color}" stroke-width="2" '
            f'stroke-dasharray="{style.dash}"/>'
        )
        endpoint_x, endpoint_y = plotted[-1]
        marks.append(
            _endpoint_label(
                sample_id,
                endpoint_x=endpoint_x,
                endpoint_y=endpoint_y,
                label_x=label_x,
                label_y=label_positions[sample_id],
                style=style,
            )
        )
        legend.append(_line_legend_entry(sample_id, assay_target, style))
    zero_reference = ""
    if x_min <= 0 <= x_max:
        zero_x = x_position(0)
        zero_reference = (
            f'<line class="zero-reference" x1="{zero_x:.1f}" y1="{plot_top:.1f}" '
            f'x2="{zero_x:.1f}" y2="{plot_bottom:.1f}" stroke="#555" '
            'stroke-width="1.5" stroke-dasharray="4 3"/>'
        )
    accessible_name = f"{title}: {x_axis_label} by {y_axis_label}"
    svg = (
        f'<svg class="panel-chart profile-chart" width="{width:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" role="img" '
        f'aria-label="{html.escape(accessible_name, quote=True)}">'
        f"<title>{html.escape(accessible_name)}</title>"
        f'<text class="axis-title axis-title-x" x="{(left + right) / 2:.1f}" '
        f'y="{height - 8:.1f}" text-anchor="middle">'
        f"{html.escape(x_axis_label)}</text>"
        f'<text class="axis-title axis-title-y" x="16" '
        f'y="{(plot_top + plot_bottom) / 2:.1f}" '
        f'transform="rotate(-90 16 {(plot_top + plot_bottom) / 2:.1f})" '
        f'text-anchor="middle">{html.escape(y_axis_label)}</text>'
        + _svg_axis_x_ticks(
            x_ticks, left=left, right=right, plot_bottom=plot_bottom,
            domain_min=x_min, domain_max=x_max,
        )
        + _svg_axis_y(
            y_ticks, left=left, right=right, plot_top=plot_top,
            plot_bottom=plot_bottom, domain_min=0.0, domain_max=y_max,
        )
        + zero_reference
        + "".join(marks)
        + "</svg>"
    )
    return (
        f'<article class="qc-panel"><h3>{html.escape(title)}</h3>'
        f'<div class="panel-scroll" role="region" '
        f'aria-label="{html.escape(title, quote=True)} chart" tabindex="0">'
        f"{svg}</div><ol class=\"series-legend compact\">"
        + "".join(legend)
        + "</ol></article>"
    )


def render_panel_grid(panels: Sequence[str], *, aria_label: str) -> str:
    """Group dashboard panels into a responsive, accessibly named grid."""
    return (
        f'<div class="panel-grid" role="group" '
        f'aria-label="{html.escape(aria_label, quote=True)}">'
        + "".join(panels)
        + "</div>"
    )
