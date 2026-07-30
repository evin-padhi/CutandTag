#!/usr/bin/env python3
"""Render a deterministic representative dashboard for visual QA."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import qc_dashboard as qc


SAMPLE_GROUPS = (
    ("701", "701_IgG", "701_CTCF", "704_GATA1", "704_RUNX1"),
    ("702", "702_IgG", "702_CTCF", "705_GATA1", "705_RUNX1"),
    ("703", "703_IgG", "703_CTCF", "706_GATA1", "706_RUNX1"),
)
TARGET_BY_SAMPLE = {
    sample_id: target
    for target, sample_ids in (
        ("IgG", ("701_IgG", "702_IgG", "703_IgG")),
        ("CTCF", ("701_CTCF", "702_CTCF", "703_CTCF")),
        ("GATA1", ("704_GATA1", "705_GATA1", "706_GATA1")),
        ("RUNX1", ("704_RUNX1", "705_RUNX1", "706_RUNX1")),
    )
    for sample_id in sample_ids
}


def _metadata() -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    for group_id, control_id, ctcf_id, gata1_id, runx1_id in SAMPLE_GROUPS:
        input_group = f"input_{group_id}"
        records[control_id] = {
            "sample_id": control_id,
            "library_id": f"L{group_id}",
            "input_group": input_group,
            "assay_target": "IgG",
            "is_control": True,
            "control_id": None,
            "expected_motif": None,
        }
        for sample_id, target, library_id in (
            (ctcf_id, "CTCF", f"L{group_id}"),
            (gata1_id, "GATA1", f"L{int(group_id) + 3}"),
            (runx1_id, "RUNX1", f"L{int(group_id) + 3}"),
        ):
            records[sample_id] = {
                "sample_id": sample_id,
                "library_id": library_id,
                "input_group": input_group,
                "assay_target": target,
                "is_control": False,
                "control_id": control_id,
                "expected_motif": target,
            }
    return records


def _demultiplex(
    metadata: dict[str, dict[str, object]],
) -> dict[str, dict[str, object]]:
    library_samples: dict[str, list[str]] = {}
    for sample_id, record in metadata.items():
        library_samples.setdefault(str(record["library_id"]), []).append(sample_id)
    records: dict[str, dict[str, object]] = {}
    for library_index, library_id in enumerate(sorted(library_samples), start=1):
        sample_ids = sorted(library_samples[library_id])
        assignments = {
            sample_id: 1_100_000 * (library_index + sample_index)
            for sample_index, sample_id in enumerate(sample_ids, start=1)
        }
        assigned = sum(assignments.values())
        ambiguous = 25_000 * library_index
        unassigned = 100_000 + 10_000 * library_index
        total = assigned + ambiguous + unassigned
        records[library_id] = {
            "total_reads": total,
            "assigned_reads": assigned,
            "ambiguous_reads": ambiguous,
            "unassigned_reads": unassigned,
            "assigned_fraction": assigned / total,
            "ambiguous_fraction": ambiguous / total,
            "unassigned_fraction": unassigned / total,
            "assignment_counts": assignments,
        }
    return records


def _libraries(
    metadata: dict[str, dict[str, object]],
) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    for index, sample_id in enumerate(sorted(metadata), start=1):
        mapped = 61.5 + 2.35 * index
        raw_reads = 2_600_000 + 735_000 * index
        mapq_fragments = 190_000 + 170_000 * index
        insert_q25 = 115 + 7 * index
        records[sample_id] = {
            "sample_id": sample_id,
            "raw_total_reads": raw_reads,
            "mapped_percent": mapped,
            "properly_paired_percent": mapped - 4.75,
            "mapq_filtered_reads": 2 * mapq_fragments,
            "mapq_filtered_fragments": mapq_fragments,
            "mapq_filtered_fraction": 0.42 + 0.025 * index,
            "markdup_examined_reads": raw_reads - 30_000,
            "duplicate_total": 105_000 + 41_000 * index,
            "duplicate_percent": 7.25 + 2.15 * index,
            "mitochondrial_percent": 0.7 + 0.3 * index,
            "estimated_library_size": 3_400_000 + 850_000 * index,
            "insert_size_total_pairs": 2_000 + 173 * index,
            "insert_size_min": 35 + index,
            "insert_size_q25": insert_q25,
            "insert_size_mean": insert_q25 + 95.25,
            "insert_size_median": insert_q25 + 68,
            "insert_size_q75": insert_q25 + 171,
            "insert_size_max": 1_050 + 23 * index,
        }
    return records


def _insert_sizes(
    metadata: dict[str, dict[str, object]],
) -> dict[str, list[dict[str, int]]]:
    boundary_sizes = (90, 249, 250, 499, 500, 749, 750, 1_001)
    return {
        sample_id: [
            {
                "insert_size": size,
                "pair_count": (index + 2) * (position + 1),
            }
            for position, size in enumerate(boundary_sizes)
        ]
        for index, sample_id in enumerate(sorted(metadata), start=1)
    }


def _peaks(
    metadata: dict[str, dict[str, object]],
) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    target_ids = [
        sample_id
        for sample_id in sorted(metadata)
        if TARGET_BY_SAMPLE[sample_id] != "IgG"
    ]
    for index, sample_id in enumerate(target_ids, start=1):
        peak_count = 8_500 + 2_850 * index
        minimum = 115 + 9 * index
        q25 = 225 + 13 * index
        median = 410 + 47 * index
        q75 = 790 + 83 * index
        maximum = 1_700 + 215 * index
        total_fragments = 420_000 + 530_000 * index
        frip = 0.12 + 0.035 * index
        records[sample_id] = {
            "peak_count": peak_count,
            "total_covered_bases": 7_500_000 + 4_375_000 * index,
            "width_count": peak_count,
            "width_min": minimum,
            "width_q25": q25,
            "width_mean": median + 115.5,
            "width_median": median,
            "width_q75": q75,
            "width_max": maximum,
            "total_fragments": total_fragments,
            "fragments_in_peaks": round(total_fragments * frip),
            "frip": frip,
        }
    return records


def _peak_widths(
    peaks: dict[str, dict[str, object]],
) -> dict[str, list[dict[str, int]]]:
    boundary_widths = (125, 249, 250, 499, 500, 749, 750, 1_250)
    return {
        sample_id: [
            {
                "width": width,
                "peak_count": (index + 1) * (len(boundary_widths) - position),
            }
            for position, width in enumerate(boundary_widths)
        ]
        for index, sample_id in enumerate(sorted(peaks), start=1)
    }


def _fragments_per_peak(
    peaks: dict[str, dict[str, object]],
) -> dict[str, list[dict[str, int]]]:
    return {
        sample_id: [
            {"fragment_count": 0, "peak_count": index + 2},
            {"fragment_count": 1, "peak_count": 34 + 2 * index},
            {"fragment_count": 2, "peak_count": 24 + 3 * index},
            {"fragment_count": 7, "peak_count": 12 + index},
            {"fragment_count": 25, "peak_count": 5 + index},
            {"fragment_count": 120, "peak_count": index},
        ]
        for index, sample_id in enumerate(sorted(peaks), start=1)
    }


def _tss(
    metadata: dict[str, dict[str, object]],
) -> dict[str, dict[str, object]]:
    positions = (-3_000, -2_000, -1_000, -250, 0, 250, 1_000, 2_000, 3_000)
    records: dict[str, dict[str, object]] = {}
    for index, sample_id in enumerate(sorted(metadata), start=1):
        target = TARGET_BY_SAMPLE[sample_id]
        target_peak = {"IgG": 2.2, "CTCF": 7.4, "GATA1": 4.8, "RUNX1": 8.1}[target]
        profile = []
        for position in positions:
            center_weight = max(0.0, 1.0 - abs(position) / 1_250)
            signal = 1.0 + target_peak * center_weight + 0.035 * index
            if abs(position) == 3_000:
                signal = 1.0 + 0.005 * (index % 3)
            profile.append((position, signal))
        records[sample_id] = {
            "status": "computed",
            "enrichment": 1.75 + 0.43 * index,
            "profile": profile,
            "score_status": "computed",
        }
    return records


def _motif_records(
    metadata: dict[str, dict[str, object]],
) -> tuple[
    dict[str, dict[str, object]],
    dict[str, list[dict[str, object]]],
    dict[str, list[dict[str, object]]],
]:
    motifs: dict[str, dict[str, object]] = {}
    top_motifs: dict[str, list[dict[str, object]]] = {}
    ame_motifs: dict[str, list[dict[str, object]]] = {}
    motif_catalog = {
        "CTCF": (
            ("MA0139.1", "CTCF", 0.0),
            ("MA1102.1", "CTCF::CTCFL", 2.5e-8),
            ("MA1929.1", "CTCFL", 0.24),
        ),
        "GATA1": (
            ("MA0035.4", "GATA1", 4.0e-12),
            ("MA0140.3", "GATA1::TAL1", 7.5e-7),
            ("MA1631.1", "TAL1::GATA1", 0.18),
        ),
        "RUNX1": (
            ("MA0002.3", "RUNX1", 9.0e-10),
            ("MA0511.2", "RUNX1::CBFB", 1.0e-5),
            ("MA0684.2", "RUNX2", 0.31),
        ),
    }
    background = (
        ("MA0098.3", "ETS1", 3.0e-18),
        ("MA0473.4", "ELF1", 2.0e-6),
        ("MA0079.5", "SP1", 0.075),
        ("MA0491.3", "JUND", 0.42),
    )
    for sample_index, sample_id in enumerate(sorted(metadata), start=1):
        target = TARGET_BY_SAMPLE[sample_id]
        if target == "IgG":
            continue
        rows = []
        for rank, (motif_id, alternate, adjusted) in enumerate(
            (*motif_catalog[target], *background),
            start=1,
        ):
            scaled_adjusted = (
                adjusted
                if adjusted in {0.0, 0.18, 0.24, 0.31, 0.42}
                else min(1.0, adjusted * (1 + sample_index / 20))
            )
            rows.append({
                "rank": rank,
                "motif_id": motif_id,
                "motif_alt_id": alternate,
                "adjusted_p_value": scaled_adjusted,
                "p_value": scaled_adjusted / 2 if scaled_adjusted else 0.0,
                "effect": 1.15 + 0.07 * rank + 0.01 * sample_index,
                "positive_sequences": 85 + 3 * rank + sample_index,
            })
        ame_motifs[sample_id] = rows
        top_motifs[sample_id] = rows[: qc.TOP_MOTIF_LIMIT]
        best = rows[0]
        motifs[sample_id] = {
            "status": "computed",
            "best_motif_id": best["motif_id"],
            "best_adjusted_p_value": best["adjusted_p_value"],
            "ame_status": "computed",
        }
    return motifs, top_motifs, ame_motifs


def build_fixture_data() -> dict[str, object]:
    """Build the representative joined report through the public interface."""
    metadata = _metadata()
    peaks = _peaks(metadata)
    motifs, top_motifs, ame_motifs = _motif_records(metadata)
    return qc.build_report_data(
        metadata,
        _demultiplex(metadata),
        _libraries(metadata),
        peaks,
        _tss(metadata),
        motifs,
        top_motifs,
        annotation_status="computed_gtf",
        insert_sizes=_insert_sizes(metadata),
        peak_widths=_peak_widths(peaks),
        fragments_per_peak=_fragments_per_peak(peaks),
        ame_motifs=ame_motifs,
        motif_analysis_status="computed",
    )


def render_fixture(output_path: Path) -> None:
    """Render the representative fixture with the public dashboard renderer."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        qc.render_dashboard(build_fixture_data()),
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Render the representative consolidated QC dashboard fixture."
    )
    parser.add_argument("output", type=Path, help="output HTML path")
    arguments = parser.parse_args(argv)
    render_fixture(arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
