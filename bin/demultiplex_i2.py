#!/usr/bin/env python3
"""Split synchronized R1/R2 FASTQs by their I2 barcode."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import tempfile
from collections import Counter
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO


class DemultiplexError(ValueError):
    """Raised when FASTQ inputs or demultiplexing parameters are invalid."""


def hamming(a: str, b: str) -> int:
    """Return the Hamming distance between two equal-length strings."""

    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(left != right for left, right in zip(a, b, strict=True))


def assign_barcode(
    observed: str, expected: dict[str, str], max_mismatches: int
) -> tuple[str, str | None]:
    """Assign an observed barcode to one unique closest expected barcode."""

    if max_mismatches < 0:
        raise ValueError("max_mismatches must be non-negative")
    candidates: list[tuple[int, str]] = []
    for sample_id, barcode in expected.items():
        if len(observed) == len(barcode):
            candidates.append((hamming(observed, barcode), sample_id))
    if not candidates:
        return "unassigned", None

    closest_distance = min(distance for distance, _ in candidates)
    if closest_distance > max_mismatches:
        return "unassigned", None
    closest_samples = [
        sample_id for distance, sample_id in candidates if distance == closest_distance
    ]
    if len(closest_samples) != 1:
        return "ambiguous", None
    return "assigned", closest_samples[0]


def _open_fastq(path: Path) -> TextIO:
    with path.open("rb") as probe:
        gzip_magic = probe.read(2) == b"\x1f\x8b"
    if gzip_magic:
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def _read_fastq_records(path: Path, label: str) -> Iterator[tuple[str, str, str, str]]:
    with _open_fastq(path) as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if header == "":
                return
            record_number += 1
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if "" in (sequence, separator, quality):
                raise DemultiplexError(
                    f"{label}: truncated FASTQ at record {record_number}"
                )
            if not header.startswith("@"):
                raise DemultiplexError(
                    f"{label}: invalid FASTQ header at record {record_number}"
                )
            if not separator.startswith("+"):
                raise DemultiplexError(
                    f"{label}: invalid FASTQ separator at record {record_number}"
                )
            sequence_value = sequence.rstrip("\r\n")
            quality_value = quality.rstrip("\r\n")
            if len(sequence_value) != len(quality_value):
                raise DemultiplexError(
                    f"{label}: sequence and quality lengths differ at record {record_number}"
                )
            yield header, sequence, separator, quality


def _normalized_read_id(header: str) -> str:
    identifier = header.rstrip("\r\n").split(maxsplit=1)[0]
    if identifier.startswith("@"):
        identifier = identifier[1:]
    return re.sub(r"/\d+$", "", identifier)


def _parse_samples(values: list[str]) -> dict[str, str]:
    expected: dict[str, str] = {}
    for value in values:
        sample_id, separator, barcode = value.partition("=")
        sample_id = sample_id.strip()
        barcode = barcode.strip().upper()
        if not separator or not sample_id or not barcode:
            raise DemultiplexError("--sample must use SAMPLE_ID=BARCODE")
        if sample_id in expected:
            raise DemultiplexError(f"duplicate sample ID {sample_id!r}")
        if barcode in expected.values():
            raise DemultiplexError(f"duplicate expected barcode {barcode!r}")
        expected[sample_id] = barcode
    if not expected:
        raise DemultiplexError("at least one --sample is required")
    barcode_lengths = {len(barcode) for barcode in expected.values()}
    if len(barcode_lengths) != 1:
        raise DemultiplexError("all expected barcodes must have equal length")
    return expected


def _atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
        Path(temporary_name).replace(path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def _metrics_paths(metrics: Path) -> tuple[Path, Path]:
    json_path = metrics if metrics.suffix == ".json" else metrics.with_suffix(".json")
    return json_path, json_path.with_suffix(".tsv")


def _write_metrics(
    metrics: Path,
    total_reads: int,
    assignment_counts: Counter[str],
    ambiguous_reads: int,
    unassigned_reads: int,
    observed_i2_counts: Counter[str],
) -> None:
    assigned_reads = sum(assignment_counts.values())
    denominator = total_reads or 1
    payload = {
        "total_reads": total_reads,
        "assigned_reads": assigned_reads,
        "ambiguous_reads": ambiguous_reads,
        "unassigned_reads": unassigned_reads,
        "assigned_fraction": assigned_reads / denominator,
        "ambiguous_fraction": ambiguous_reads / denominator,
        "unassigned_fraction": unassigned_reads / denominator,
        "assignment_counts": dict(sorted(assignment_counts.items())),
        "observed_i2_counts": dict(sorted(observed_i2_counts.items())),
    }
    lines = ["metric\tvalue"]
    for name in (
        "total_reads",
        "assigned_reads",
        "ambiguous_reads",
        "unassigned_reads",
        "assigned_fraction",
        "ambiguous_fraction",
        "unassigned_fraction",
    ):
        lines.append(f"{name}\t{payload[name]}")
    for sample_id, count in payload["assignment_counts"].items():
        lines.append(f"assigned:{sample_id}\t{count}")
    for barcode, count in payload["observed_i2_counts"].items():
        lines.append(f"observed_i2:{barcode}\t{count}")
    json_path, tsv_path = _metrics_paths(metrics)
    _atomic_write(json_path, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    _atomic_write(tsv_path, "\n".join(lines) + "\n")


def _temporary_backup_path(path: Path) -> Path:
    descriptor, name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".backup", dir=path.parent
    )
    os.close(descriptor)
    backup = Path(name)
    backup.unlink()
    return backup


def _publish_staged_outputs(staged_outputs: list[tuple[Path, Path]]) -> None:
    """Atomically publish staged files, rolling back already-replaced paths on error."""

    stage_device = os.stat(staged_outputs[0][0].parent).st_dev
    for _, final_path in staged_outputs:
        final_path.parent.mkdir(parents=True, exist_ok=True)
        if os.stat(final_path.parent).st_dev != stage_device:
            raise DemultiplexError(
                "all output paths must be on the staging filesystem for atomic publication"
            )

    published: list[tuple[Path, Path | None]] = []
    try:
        for staged_path, final_path in staged_outputs:
            backup_path = None
            if final_path.exists():
                backup_path = _temporary_backup_path(final_path)
                os.replace(final_path, backup_path)
            published.append((final_path, backup_path))
            os.replace(staged_path, final_path)
    except BaseException:
        for final_path, backup_path in reversed(published):
            final_path.unlink(missing_ok=True)
            if backup_path is not None and backup_path.exists():
                os.replace(backup_path, final_path)
        raise
    else:
        for _, backup_path in published:
            if backup_path is not None:
                try:
                    backup_path.unlink(missing_ok=True)
                except OSError:
                    pass


def demultiplex(
    r1_path: Path,
    r2_path: Path,
    i2_path: Path,
    expected: dict[str, str],
    max_mismatches: int,
    outdir: Path,
    metrics: Path,
    allow_empty: bool,
) -> None:
    """Stream and split three synchronized FASTQs, then emit output metrics."""

    barcode_length = len(next(iter(expected.values())))
    if not 0 <= max_mismatches < barcode_length:
        raise DemultiplexError(
            "max mismatches must be non-negative and smaller than barcode length"
        )
    for path, label in ((r1_path, "R1"), (r2_path, "R2"), (i2_path, "I2")):
        if not path.is_file():
            raise DemultiplexError(f"{label} FASTQ file does not exist: {path}")

    outdir.mkdir(parents=True, exist_ok=True)
    metrics_json, metrics_tsv = _metrics_paths(metrics)
    metrics_json.parent.mkdir(parents=True, exist_ok=True)
    metrics_tsv.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".demultiplex-", dir=outdir) as temporary_name:
        staging_directory = Path(temporary_name)
        staged_fastqs = staging_directory / "fastqs"
        staged_fastqs.mkdir()
        writers: dict[str, tuple[TextIO, TextIO]] = {}
        try:
            for sample_id in expected:
                writers[sample_id] = (
                    gzip.open(
                        staged_fastqs / f"{sample_id}_R1.fastq.gz", "wt", encoding="utf-8"
                    ),
                    gzip.open(
                        staged_fastqs / f"{sample_id}_R2.fastq.gz", "wt", encoding="utf-8"
                    ),
                )

            r1_records = _read_fastq_records(r1_path, "R1")
            r2_records = _read_fastq_records(r2_path, "R2")
            i2_records = _read_fastq_records(i2_path, "I2")
            sentinel = object()
            total_reads = ambiguous_reads = unassigned_reads = 0
            assignment_counts: Counter[str] = Counter({sample_id: 0 for sample_id in expected})
            observed_i2_counts: Counter[str] = Counter()

            while True:
                r1_record = next(r1_records, sentinel)
                r2_record = next(r2_records, sentinel)
                i2_record = next(i2_records, sentinel)
                if r1_record is sentinel and r2_record is sentinel and i2_record is sentinel:
                    break
                record_number = total_reads + 1
                if sentinel in (r1_record, r2_record, i2_record):
                    raise DemultiplexError(
                        f"FASTQ record counts differ at record {record_number}"
                    )
                assert isinstance(r1_record, tuple)
                assert isinstance(r2_record, tuple)
                assert isinstance(i2_record, tuple)
                identifiers = tuple(
                    _normalized_read_id(record[0])
                    for record in (r1_record, r2_record, i2_record)
                )
                if len(set(identifiers)) != 1:
                    raise DemultiplexError(
                        f"read identifiers differ at record {record_number}: {identifiers}"
                    )
                observed = i2_record[1].rstrip("\r\n").upper()
                observed_i2_counts[observed] += 1
                status, sample_id = assign_barcode(observed, expected, max_mismatches)
                total_reads += 1
                if status == "assigned":
                    assert sample_id is not None
                    r1_writer, r2_writer = writers[sample_id]
                    r1_writer.write("".join(r1_record))
                    r2_writer.write("".join(r2_record))
                    assignment_counts[sample_id] += 1
                elif status == "ambiguous":
                    ambiguous_reads += 1
                else:
                    unassigned_reads += 1
        finally:
            for r1_writer, r2_writer in writers.values():
                r1_writer.close()
                r2_writer.close()

        empty_samples = [
            sample_id for sample_id, count in assignment_counts.items() if count == 0
        ]
        if empty_samples and not allow_empty:
            raise DemultiplexError(
                f"empty derived samples: {', '.join(empty_samples)}; use --allow-empty to permit"
            )
        staged_metrics = staging_directory / "metrics.json"
        _write_metrics(
            staged_metrics,
            total_reads,
            assignment_counts,
            ambiguous_reads,
            unassigned_reads,
            observed_i2_counts,
        )
        staged_outputs = [
            (staged_fastqs / f"{sample_id}_{read}.fastq.gz", outdir / f"{sample_id}_{read}.fastq.gz")
            for sample_id in expected
            for read in ("R1", "R2")
        ]
        staged_outputs.extend(
            [
                (staged_metrics, metrics_json),
                (staged_metrics.with_suffix(".tsv"), metrics_tsv),
            ]
        )
        _publish_staged_outputs(staged_outputs)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--i2", type=Path, required=True)
    parser.add_argument("--sample", action="append", required=True, metavar="SAMPLE_ID=BARCODE")
    parser.add_argument("--max-mismatches", type=int, default=0)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--allow-empty", action="store_true")
    arguments = parser.parse_args(argv)
    try:
        expected = _parse_samples(arguments.sample)
        demultiplex(
            arguments.r1,
            arguments.r2,
            arguments.i2,
            expected,
            arguments.max_mismatches,
            arguments.outdir,
            arguments.metrics,
            arguments.allow_empty,
        )
    except DemultiplexError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
