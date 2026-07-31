#!/usr/bin/env python3
"""Interval parsing and peak-enrichment sampling/statistics helpers."""

from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
import math
from pathlib import Path
import random
from statistics import mean, pstdev
from typing import Iterable, Sequence


@dataclass(frozen=True, slots=True)
class Interval:
    chrom: str
    start: int
    end: int

    @property
    def width(self) -> int:
        return self.end - self.start


@dataclass(frozen=True, slots=True)
class ReferenceRecord:
    reference_id: str
    tf: str
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
) -> list[Interval] | None:
    foreground = list(foreground)
    blacklist = list(blacklist)
    if not foreground:
        return []
    if max_attempts <= 0:
        return None

    placed: list[Interval] = []
    gc_models = {"gc_matched", "length_gc_matched"}

    for foreground_interval in foreground:
        target_gc = gc_fraction(foreground_interval, fasta_sequences) if model in gc_models else None
        attempts_1pct = min(max_attempts, max(1, max_attempts // 2)) if model in gc_models else max_attempts
        attempts_fallback = max_attempts - attempts_1pct if model in gc_models else 0

        found: Interval | None = None
        for tolerance, attempts in ((0.01, attempts_1pct), (None if target_gc is None else 0.02, attempts_fallback)):
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
    seed: int,
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
