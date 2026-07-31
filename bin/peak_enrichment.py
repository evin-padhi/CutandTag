#!/usr/bin/env python3
"""Interval parsing and peak-enrichment sampling/statistics helpers."""

from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import defaultdict
import csv
from dataclasses import dataclass
from pathlib import Path
import random
import re
from statistics import mean, pstdev
import sys
from typing import Iterable, Sequence


IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
ENRICHMENT_FIELDNAMES = [
    "foreground_id",
    "foreground_tf",
    "reference_id",
    "reference_tf",
    "background_model",
    "seed",
    "status",
    "permutations_requested",
    "permutations_succeeded",
    "observed_overlap_count",
    "null_mean_overlap_count",
    "null_stddev_overlap_count",
    "enrichment_ratio",
    "empirical_p_value",
]
STATUS_FIELDNAMES = [
    "foreground_id",
    "foreground_tf",
    "reference_id",
    "reference_tf",
    "background_model",
    "seed",
    "status",
    "permutations_requested",
    "permutations_succeeded",
]


@dataclass(frozen=True, slots=True)
class Interval:
    chrom: str
    start: int
    end: int

    @property
    def width(self) -> int:
        return self.end - self.start


@dataclass(frozen=True, slots=True)
class ForegroundRecord:
    foreground_id: str
    foreground_tf: str
    peak_file: Path


@dataclass(frozen=True, slots=True)
class ReferenceRecord:
    reference_id: str
    tf: str
    reference_type: str
    peak_file: Path


def _validate_interval(interval: Interval, label: str = "interval") -> None:
    if interval.start < 0 or interval.end < 0:
        raise ValueError(f"{label}: coordinates must be non-negative")
    if interval.start >= interval.end:
        raise ValueError(f"{label}: start must be less than end")


def _overlaps(left: Interval, right: Interval) -> bool:
    return left.chrom == right.chrom and left.start < right.end and left.end > right.start


def _group_by_chrom(intervals: Iterable[Interval]) -> dict[str, list[Interval]]:
    grouped: dict[str, list[Interval]] = defaultdict(list)
    for interval in intervals:
        grouped[interval.chrom].append(interval)
    return grouped


def parse_intervals(path: str | Path, chrom_sizes: dict[str, int]) -> list[Interval]:
    intervals: list[Interval] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"line {line_number}: expected at least three BED columns")
            chrom = fields[0]
            if chrom not in chrom_sizes:
                raise ValueError(f"line {line_number}: unknown chromosome {chrom!r}")
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as error:
                raise ValueError(f"line {line_number}: coordinates must be integers") from error
            interval = Interval(chrom, start, end)
            _validate_interval(interval, f"line {line_number}")
            if end > chrom_sizes[chrom]:
                raise ValueError(f"line {line_number}: end exceeds chromosome size")
            intervals.append(interval)
    return intervals


def count_overlapping_foreground(
    foreground: Sequence[Interval],
    reference: Sequence[Interval],
) -> int:
    foreground_by_chrom = _group_by_chrom(foreground)
    reference_by_chrom = _group_by_chrom(reference)
    overlap_count = 0

    for chrom, foreground_intervals in foreground_by_chrom.items():
        references = sorted(reference_by_chrom.get(chrom, ()), key=lambda interval: (interval.start, interval.end))
        if not references:
            continue
        starts = [interval.start for interval in references]
        prefix_max_ends: list[int] = []
        running_max = -1
        for interval in references:
            running_max = max(running_max, interval.end)
            prefix_max_ends.append(running_max)
        for interval in sorted(foreground_intervals, key=lambda item: (item.start, item.end)):
            hi = bisect_left(starts, interval.end)
            if hi and prefix_max_ends[hi - 1] > interval.start:
                overlap_count += 1

    return overlap_count


def gc_fraction(interval: Interval, fasta_sequences: dict[str, str]) -> float:
    _validate_interval(interval)
    try:
        sequence = fasta_sequences[interval.chrom]
    except KeyError as error:
        raise ValueError(f"unknown chromosome {interval.chrom!r}") from error
    if interval.end > len(sequence):
        raise ValueError(f"{interval.chrom}: interval extends beyond the available sequence")
    window = sequence[interval.start : interval.end]
    gc_bases = sum(base.upper() in {"G", "C"} for base in window)
    return gc_bases / len(window)


def _candidate_overlaps(candidate: Interval, others: Sequence[Interval]) -> bool:
    return any(_overlaps(candidate, other) for other in others if other.chrom == candidate.chrom)


def _pick_candidate_coordinates(
    foreground: Sequence[Interval],
    chrom_sizes: dict[str, int],
    rng: random.Random,
    model: str,
    foreground_interval: Interval,
) -> tuple[str, int]:
    foreground_chroms = [interval.chrom for interval in foreground]
    foreground_widths = [interval.width for interval in foreground]

    if model in {"length_matched", "length_gc_matched"}:
        chrom = foreground_interval.chrom
        width = foreground_interval.width
    elif model in {"random", "gc_matched"}:
        chrom = rng.choice(foreground_chroms)
        width = rng.choice(foreground_widths)
    else:
        raise ValueError(f"unknown background model {model!r}")

    if chrom not in chrom_sizes:
        raise ValueError(f"chromosome {chrom!r} is missing from the chromosome size map")
    return chrom, width


def _sample_candidate(
    foreground: Sequence[Interval],
    chrom_sizes: dict[str, int],
    rng: random.Random,
    model: str,
    foreground_interval: Interval,
) -> Interval | None:
    chrom, width = _pick_candidate_coordinates(foreground, chrom_sizes, rng, model, foreground_interval)
    if width > chrom_sizes[chrom]:
        return None
    start = rng.randrange(0, chrom_sizes[chrom] - width + 1)
    return Interval(chrom, start, start + width)


def sample_background(
    foreground: Sequence[Interval],
    chrom_sizes: dict[str, int],
    fasta_sequences: dict[str, str],
    blacklist: Sequence[Interval],
    model: str,
    rng: random.Random,
    max_attempts: int,
    gc_tolerance: float = 0.02,
) -> list[Interval] | None:
    foreground = list(foreground)
    blacklist = list(blacklist)
    if not foreground:
        return []
    if max_attempts <= 0:
        return None
    if gc_tolerance <= 0:
        raise ValueError("gc_tolerance must be positive")

    placed: list[Interval] = []
    gc_models = {"gc_matched", "length_gc_matched"}

    for foreground_interval in foreground:
        target_gc = gc_fraction(foreground_interval, fasta_sequences) if model in gc_models else None
        attempts_1pct = min(max_attempts, max(1, max_attempts // 2)) if model in gc_models else max_attempts
        attempts_fallback = max_attempts - attempts_1pct if model in gc_models else 0

        found: Interval | None = None
        for tolerance, attempts in ((0.01, attempts_1pct), (None if target_gc is None else gc_tolerance, attempts_fallback)):
            for _ in range(attempts):
                candidate = _sample_candidate(foreground, chrom_sizes, rng, model, foreground_interval)
                if candidate is None:
                    continue
                if _candidate_overlaps(candidate, blacklist):
                    continue
                if _candidate_overlaps(candidate, placed):
                    continue
                if target_gc is not None:
                    candidate_gc = gc_fraction(candidate, fasta_sequences)
                    if abs(candidate_gc - target_gc) > tolerance:
                        continue
                found = candidate
                break
            if found is not None:
                break

        if found is None:
            return None
        placed.append(found)

    return placed


def _blank_row(model: str, seed: str, status: str, permutations: int) -> dict[str, object]:
    return {
        "background_model": model,
        "seed": seed,
        "status": status,
        "permutations_requested": permutations,
        "permutations_succeeded": 0,
        "observed_overlap_count": None,
        "null_mean_overlap_count": None,
        "null_stddev_overlap_count": None,
        "enrichment_ratio": None,
        "empirical_p_value": None,
    }


def calculate_enrichment(
    foreground: Sequence[Interval],
    reference: Sequence[Interval],
    chrom_sizes: dict[str, int],
    fasta_sequences: dict[str, str],
    blacklist: Sequence[Interval],
    permutations: int,
    seed: int | str,
    gc_tolerance: float,
) -> list[dict[str, object]]:
    if permutations < 0:
        raise ValueError("permutations must be non-negative")
    if gc_tolerance <= 0:
        raise ValueError("gc_tolerance must be positive")

    foreground = list(foreground)
    reference = list(reference)
    blacklist = list(blacklist)
    models = ("random", "length_matched", "gc_matched", "length_gc_matched")
    rows: list[dict[str, object]] = []

    for model in models:
        model_seed = f"{seed}:{model}"
        if not foreground:
            rows.append(_blank_row(model, model_seed, "no_foreground_peaks", permutations))
            continue
        if not reference:
            rows.append(_blank_row(model, model_seed, "no_reference_peaks", permutations))
            continue

        model_rng = random.Random(model_seed)
        observed_overlap_count = count_overlapping_foreground(foreground, reference)
        null_counts: list[int] = []
        for _ in range(permutations):
            background = sample_background(
                foreground,
                chrom_sizes,
                fasta_sequences,
                blacklist,
                model,
                model_rng,
                max_attempts=1000,
                gc_tolerance=gc_tolerance,
            )
            if background is None:
                continue
            null_counts.append(count_overlapping_foreground(background, reference))

        successful = len(null_counts)
        if successful == 0:
            rows.append(_blank_row(model, model_seed, "insufficient_background", permutations))
            rows[-1]["observed_overlap_count"] = observed_overlap_count
            continue

        null_mean_overlap_count = mean(null_counts)
        null_stddev_overlap_count = pstdev(null_counts) if successful > 1 else 0.0
        extreme_nulls = sum(count >= observed_overlap_count for count in null_counts)
        empirical_p_value = (1 + extreme_nulls) / (1 + successful)
        enrichment_ratio = (
            observed_overlap_count / null_mean_overlap_count if null_mean_overlap_count else None
        )

        rows.append(
            {
                "background_model": model,
                "seed": model_seed,
                "status": "ok",
                "permutations_requested": permutations,
                "permutations_succeeded": successful,
                "observed_overlap_count": observed_overlap_count,
                "null_mean_overlap_count": null_mean_overlap_count,
                "null_stddev_overlap_count": null_stddev_overlap_count,
                "enrichment_ratio": enrichment_ratio,
                "empirical_p_value": empirical_p_value,
            }
        )

    return rows


def load_fasta_sequences(path: str | Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current_name: str | None = None
    fasta_path = Path(path)
    if not fasta_path.is_file():
        raise ValueError(f"FASTA file does not exist: {fasta_path}")
    for raw_line in fasta_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current_name = line[1:].split()[0]
            if not current_name:
                raise ValueError("FASTA record is missing a sequence name")
            sequences[current_name] = []
            continue
        if current_name is None:
            raise ValueError("FASTA sequence data appeared before the first header")
        sequences[current_name].append(line)
    if not sequences:
        raise ValueError("FASTA contains no records")
    return {name: "".join(parts) for name, parts in sequences.items()}


def _manifest_dialect(sample: str) -> csv.Dialect:
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t")
    except csv.Error:
        return csv.excel_tab if "\t" in sample.partition("\n")[0] else csv.excel


def _read_manifest_rows(
    path: str | Path,
    required_columns: Sequence[str],
    optional_columns: Sequence[str] = (),
    allow_empty: bool = False,
) -> list[tuple[int, dict[str, str]]]:
    manifest_path = Path(path).resolve()
    if not manifest_path.is_file():
        raise ValueError(f"manifest file does not exist: {manifest_path}")

    with manifest_path.open(newline="", encoding="utf-8") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        reader = csv.DictReader(handle, dialect=_manifest_dialect(sample))
        fieldnames = reader.fieldnames or []
        missing_columns = [column for column in required_columns if column not in fieldnames]
        if missing_columns:
            missing = ", ".join(missing_columns)
            raise ValueError(f"manifest is missing required columns: {missing}")

        rows: list[tuple[int, dict[str, str]]] = []
        for row_number, row in enumerate(reader, start=2):
            normalized: dict[str, str] = {}
            for column in required_columns:
                value = (row.get(column) or "").strip()
                if not value:
                    raise ValueError(f"row {row_number}: blank required value for {column}")
                normalized[column] = value
            for column in optional_columns:
                if column in fieldnames:
                    normalized[column] = (row.get(column) or "").strip()
            rows.append((row_number, normalized))

    if not rows and not allow_empty:
        raise ValueError("manifest contains no rows")
    return rows


def _validate_identifier(value: str, label: str, row_number: int) -> None:
    if IDENTIFIER_PATTERN.fullmatch(value) is None:
        raise ValueError(
            f"row {row_number}: invalid {label} {value!r}; "
            "use only letters, digits, dot, underscore, or dash, "
            "and begin with a letter or digit"
        )


def _resolve_manifest_path(value: str, manifest_path: Path, row_number: int, column: str) -> Path:
    candidate = Path(value)
    if not candidate.is_absolute():
        candidate = manifest_path.parent / candidate
    resolved = candidate.resolve()
    if not resolved.is_file():
        raise ValueError(f"row {row_number}: {column} file does not exist: {resolved}")
    return resolved


def load_foreground_manifest(path: str | Path, chrom_sizes: dict[str, int]) -> list[ForegroundRecord]:
    manifest_path = Path(path).resolve()
    rows = _read_manifest_rows(
        manifest_path,
        ("foreground_id", "foreground_tf", "peak_file"),
        allow_empty=True,
    )
    records: list[ForegroundRecord] = []
    seen_ids: set[str] = set()

    for row_number, row in rows:
        foreground_id = row["foreground_id"]
        _validate_identifier(foreground_id, "foreground_id", row_number)
        if foreground_id in seen_ids:
            raise ValueError(f"row {row_number}: duplicate foreground_id {foreground_id!r}")
        seen_ids.add(foreground_id)
        peak_file = _resolve_manifest_path(row["peak_file"], manifest_path, row_number, "peak_file")
        parse_intervals(peak_file, chrom_sizes)
        records.append(ForegroundRecord(foreground_id, row["foreground_tf"], peak_file))

    return records


def load_reference_manifest(path: str | Path, chrom_sizes: dict[str, int]) -> list[ReferenceRecord]:
    manifest_path = Path(path).resolve()
    rows = _read_manifest_rows(
        manifest_path,
        ("reference_id", "tf", "peak_file"),
        optional_columns=("reference_type",),
    )
    records: list[ReferenceRecord] = []
    seen_ids: set[str] = set()

    for row_number, row in rows:
        reference_id = row["reference_id"]
        _validate_identifier(reference_id, "reference_id", row_number)
        if reference_id in seen_ids:
            raise ValueError(f"row {row_number}: duplicate reference_id {reference_id!r}")
        seen_ids.add(reference_id)
        reference_type = row.get("reference_type") or "chipseq"
        if reference_type not in {"chipseq", "called_tf"}:
            raise ValueError(
                f"row {row_number}: unsupported reference_type {reference_type!r}"
            )
        peak_file = _resolve_manifest_path(row["peak_file"], manifest_path, row_number, "peak_file")
        parse_intervals(peak_file, chrom_sizes)
        records.append(ReferenceRecord(reference_id, row["tf"], reference_type, peak_file))

    return records


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _stringify(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def _write_matrix_files(rows: Sequence[dict[str, object]], outdir: Path) -> None:
    models = ("random", "length_matched", "gc_matched", "length_gc_matched")
    foreground_ids = sorted({str(row["foreground_id"]) for row in rows})
    reference_ids = sorted({str(row["reference_id"]) for row in rows})

    for model in models:
        model_rows = [row for row in rows if row["background_model"] == model]
        ratio_lookup = {
            (str(row["foreground_id"]), str(row["reference_id"])): row["enrichment_ratio"]
            for row in model_rows
            if row["status"] == "ok"
        }
        matrix_records: list[dict[str, object]] = []
        heatmap_values: list[list[float]] = []
        for foreground_id in foreground_ids:
            matrix_row: dict[str, object] = {"foreground_id": foreground_id}
            numeric_row: list[float] = []
            for reference_id in reference_ids:
                value = ratio_lookup.get((foreground_id, reference_id))
                matrix_row[reference_id] = _stringify(value)
                numeric_row.append(float("nan") if value is None else float(value))
            matrix_records.append(matrix_row)
            heatmap_values.append(numeric_row)

        _write_tsv(
            outdir / f"matrix_{model}.tsv",
            ["foreground_id", *reference_ids],
            matrix_records,
        )
        _write_heatmap_png(
            outdir / f"matrix_{model}.png",
            heatmap_values,
            foreground_ids,
            reference_ids,
            title=f"Peak enrichment ({model})",
        )


def _load_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot

    return pyplot


def _write_heatmap_png(
    path: Path,
    values: Sequence[Sequence[float]],
    row_labels: Sequence[str],
    column_labels: Sequence[str],
    title: str,
) -> None:
    pyplot = _load_pyplot()
    figure, axis = pyplot.subplots(
        figsize=(max(4, len(column_labels) * 1.2), max(3, len(row_labels) * 0.8))
    )
    if row_labels and column_labels:
        image = axis.imshow(values, aspect="auto", cmap="viridis")
        axis.set_xticks(range(len(column_labels)))
        axis.set_xticklabels(column_labels, rotation=45, ha="right")
        axis.set_yticks(range(len(row_labels)))
        axis.set_yticklabels(row_labels)
        figure.colorbar(image, ax=axis, label="Enrichment ratio")
    else:
        axis.text(0.5, 0.5, "No comparisons available", ha="center", va="center")
        axis.set_xticks([])
        axis.set_yticks([])
    axis.set_title(title)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    pyplot.close(figure)


def _write_observed_vs_null_plot(path: Path, rows: Sequence[dict[str, object]]) -> None:
    pyplot = _load_pyplot()
    figure, axis = pyplot.subplots(figsize=(max(6, len(rows) * 1.0), 4))
    ok_rows = [row for row in rows if row["status"] == "ok"]
    if ok_rows:
        ordered_rows = sorted(
            ok_rows,
            key=lambda row: (
                str(row["foreground_id"]),
                str(row["reference_id"]),
                str(row["background_model"]),
            ),
        )
        labels = [
            f"{row['foreground_id']}→{row['reference_id']}\n{row['background_model']}"
            for row in ordered_rows
        ]
        positions = list(range(len(ordered_rows)))
        observed = [float(row["observed_overlap_count"]) for row in ordered_rows]
        null_mean = [float(row["null_mean_overlap_count"]) for row in ordered_rows]
        axis.plot(positions, observed, marker="o", label="Observed overlaps")
        axis.plot(positions, null_mean, marker="s", label="Null mean overlaps")
        axis.set_xticks(positions)
        axis.set_xticklabels(labels, rotation=45, ha="right")
        axis.legend()
    else:
        axis.text(0.5, 0.5, "No successful comparisons available", ha="center", va="center")
        axis.set_xticks([])
    axis.set_ylabel("Overlap count")
    axis.set_title("Observed vs null overlap summary")
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    pyplot.close(figure)


def _status_only_rows_for_empty_foregrounds(
    permutations: int,
    seed: int | str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for model in ("random", "length_matched", "gc_matched", "length_gc_matched"):
        rows.append(
            {
                "foreground_id": "",
                "foreground_tf": "",
                "reference_id": "",
                "reference_tf": "",
                "background_model": model,
                "seed": seed,
                "status": "no_foreground_peaks",
                "permutations_requested": permutations,
                "permutations_succeeded": 0,
            }
        )
    return rows


def run_cli(
    foreground_manifest: str | Path,
    reference_manifest: str | Path,
    fasta: str | Path,
    outdir: str | Path,
    permutations: int,
    seed: int,
    gc_tolerance: float,
    blacklist: str | Path | None = None,
) -> list[dict[str, object]]:
    if permutations < 0:
        raise ValueError("permutations must be non-negative")
    if not 0 < gc_tolerance <= 1:
        raise ValueError("gc_tolerance must be in the interval (0, 1]")

    fasta_sequences = load_fasta_sequences(fasta)
    chrom_sizes = {chrom: len(sequence) for chrom, sequence in fasta_sequences.items()}
    foreground_records = sorted(
        load_foreground_manifest(foreground_manifest, chrom_sizes),
        key=lambda record: record.foreground_id,
    )
    reference_records = sorted(
        load_reference_manifest(reference_manifest, chrom_sizes),
        key=lambda record: record.reference_id,
    )
    blacklist_intervals = parse_intervals(blacklist, chrom_sizes) if blacklist else []

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)

    if not foreground_records:
        _write_tsv(outdir_path / "peak_enrichment.tsv", ENRICHMENT_FIELDNAMES, [])
        status_rows = _status_only_rows_for_empty_foregrounds(permutations, seed)
        _write_tsv(outdir_path / "enrichment_status.tsv", STATUS_FIELDNAMES, status_rows)
        _write_matrix_files([], outdir_path)
        _write_observed_vs_null_plot(outdir_path / "observed_vs_null.png", [])
        return []

    all_rows: list[dict[str, object]] = []
    for foreground in foreground_records:
        foreground_intervals = parse_intervals(foreground.peak_file, chrom_sizes)
        for reference in reference_records:
            if (
                reference.reference_type == "called_tf"
                and reference.reference_id == foreground.foreground_id
            ):
                continue
            reference_intervals = parse_intervals(reference.peak_file, chrom_sizes)
            pair_seed = f"{seed}:{foreground.foreground_id}:{reference.reference_id}"
            enrichment_rows = calculate_enrichment(
                foreground_intervals,
                reference_intervals,
                chrom_sizes,
                fasta_sequences,
                blacklist_intervals,
                permutations=permutations,
                seed=pair_seed,
                gc_tolerance=gc_tolerance,
            )
            for row in enrichment_rows:
                all_rows.append(
                    {
                        "foreground_id": foreground.foreground_id,
                        "foreground_tf": foreground.foreground_tf,
                        "reference_id": reference.reference_id,
                        "reference_tf": reference.tf,
                        **row,
                    }
                )

    ordered_rows = sorted(
        all_rows,
        key=lambda row: (
            str(row["foreground_id"]),
            str(row["reference_id"]),
            str(row["background_model"]),
        ),
    )
    _write_tsv(outdir_path / "peak_enrichment.tsv", ENRICHMENT_FIELDNAMES, ordered_rows)
    _write_tsv(
        outdir_path / "enrichment_status.tsv",
        STATUS_FIELDNAMES,
        [{field: row[field] for field in STATUS_FIELDNAMES} for row in ordered_rows],
    )
    _write_matrix_files(ordered_rows, outdir_path)
    _write_observed_vs_null_plot(outdir_path / "observed_vs_null.png", ordered_rows)
    return ordered_rows


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--foreground-manifest", required=True)
    parser.add_argument("--reference-manifest", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--gc-tolerance", type=float, default=0.02)
    parser.add_argument("--blacklist")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_argument_parser().parse_args(argv)
    try:
        run_cli(
            foreground_manifest=args.foreground_manifest,
            reference_manifest=args.reference_manifest,
            fasta=args.fasta,
            outdir=args.outdir,
            permutations=args.permutations,
            seed=args.seed,
            gc_tolerance=args.gc_tolerance,
            blacklist=args.blacklist,
        )
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
