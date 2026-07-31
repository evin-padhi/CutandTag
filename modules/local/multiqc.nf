def demuxLibraryId(meta) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    def libraryId = meta instanceof Map ? meta.library_id?.toString() : null
    if (libraryId == null || !libraryId.matches(safeId)) {
        throw new IllegalArgumentException(
            "demultiplex library_id must match ${safeId}, got ${libraryId}"
        )
    }
    libraryId
}

def libraryQcSampleId(meta) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    def sampleId = meta instanceof Map ? meta.sample_id?.toString() : null
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "library-QC sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    sampleId
}

def motifQcSampleId(meta) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    def sampleId = meta instanceof Map ? meta.sample_id?.toString() : null
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "motif-QC sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    sampleId
}

process DEMUX_QC_CUSTOM {
    tag "${meta.library_id}"
    label 'process_light'

    conda "${projectDir}/envs/multiqc.yml"
    container 'quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0'

    input:
    tuple val(meta),
        path(demux_json, stageAs: 'demultiplex.metrics.json'),
        path(demux_tsv, stageAs: 'demultiplex.metrics.tsv')

    output:
    tuple val(meta), path("*.demultiplex_qc.tsv"), emit: custom

    script:
    def libraryId = demuxLibraryId(meta)
    def outputName = "${libraryId}.demultiplex_qc.tsv"

    """
    set -euo pipefail
    python - "${libraryId}" "${outputName}" <<'PY'
import csv
import json
import sys

library_id, output_name = sys.argv[1:]
with open("demultiplex.metrics.json", encoding="utf-8") as handle:
    metrics = json.load(handle)

columns = [
    "library_id",
    "total_read_pairs",
    "assigned_read_pairs",
    "ambiguous_read_pairs",
    "unassigned_read_pairs",
    "assigned_fraction",
    "ambiguous_fraction",
    "unassigned_fraction",
]
row = {
    "library_id": library_id,
    "total_read_pairs": metrics.get("total_reads", 0),
    "assigned_read_pairs": metrics.get("assigned_reads", 0),
    "ambiguous_read_pairs": metrics.get("ambiguous_reads", 0),
    "unassigned_read_pairs": metrics.get("unassigned_reads", 0),
    "assigned_fraction": metrics.get("assigned_fraction", 0),
    "ambiguous_fraction": metrics.get("ambiguous_fraction", 0),
    "unassigned_fraction": metrics.get("unassigned_fraction", 0),
}
with open(output_name, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\\t")
    writer.writeheader()
    writer.writerow(row)
PY
    """
}

process LIBRARY_QC_CUSTOM {
    tag "${meta.sample_id}"
    label 'process_light'

    conda "${projectDir}/envs/multiqc.yml"
    container 'quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0'

    input:
    tuple val(meta),
        path(flagstat, stageAs: 'flagstat.txt'),
        path(stats, stageAs: 'alignment.stats.txt'),
        path(idxstats, stageAs: 'idxstats.tsv'),
        path(insert_size, stageAs: 'insert_size.tsv'),
        path(duplicate_metrics, stageAs: 'duplicate_metrics.json'),
        path(filtered_bam_qc, stageAs: 'filtered_bam_qc.tsv')

    output:
    tuple val(meta),
        path("*.library_qc.tsv"),
        path("*.insert_size_distribution.tsv"),
        emit: custom

    script:
    def sampleId = libraryQcSampleId(meta)
    def outputName = "${sampleId}.library_qc.tsv"
    def insertOutputName = "${sampleId}.insert_size_distribution.tsv"

    """
    set -euo pipefail
    python - "${sampleId}" "${outputName}" "${insertOutputName}" <<'PY'
import csv
import json
import math
import re
import sys

sample_id, output_name, insert_output_name = sys.argv[1:]
flagstat = open("flagstat.txt", encoding="utf-8").read().splitlines()

def percent_for(label):
    for line in flagstat:
        if label in line:
            match = re.search(r"\\(([0-9]+(?:\\.[0-9]+)?)%", line)
            if match:
                return float(match.group(1))
    raise SystemExit(f"flagstat is missing percentage for {label.strip()}")

def require_non_negative_number(mapping, key, source="markdup"):
    if key not in mapping:
        raise SystemExit(f"missing required {source} field: {key}")
    value = mapping[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SystemExit(f"{source} field {key} must be a JSON number")
    if not math.isfinite(value) or value < 0:
        raise SystemExit(f"{source} field {key} must be finite and non-negative")
    return value

stats_values = {}
with open("alignment.stats.txt", encoding="utf-8") as handle:
    for line in handle:
        fields = line.rstrip("\\n").split("\\t")
        if len(fields) >= 3 and fields[0] == "SN":
            key = fields[1].strip().rstrip(":").lower()
            value_text = fields[2].strip()
            try:
                value = float(value_text)
            except ValueError as error:
                raise SystemExit(
                    f"alignment.stats.txt SN field {key} must be numeric"
                ) from error
            if not math.isfinite(value) or value < 0:
                raise SystemExit(
                    f"alignment.stats.txt SN field {key} must be finite and non-negative"
                )
            stats_values[key] = value
if "raw total sequences" not in stats_values:
    raise SystemExit("alignment.stats.txt is missing SN raw total sequences")
raw_total_reads = stats_values["raw total sequences"]

filtered_metrics = {}
with open("filtered_bam_qc.tsv", encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\\t"):
        filtered_metrics[row.get("metric", "")] = row.get("value", "")
for required in ("filtered_reads", "filtered_fragments"):
    try:
        value = float(filtered_metrics[required])
    except (KeyError, ValueError) as error:
        raise SystemExit(
            f"filtered_bam_qc.tsv has invalid or missing {required}"
        ) from error
    if not math.isfinite(value) or value < 0 or not value.is_integer():
        raise SystemExit(
            f"filtered_bam_qc.tsv {required} must be a non-negative integer"
        )
    filtered_metrics[required] = int(value)
filtered_reads = filtered_metrics["filtered_reads"]
filtered_fragments = filtered_metrics["filtered_fragments"]
if filtered_reads != filtered_fragments * 2:
    raise SystemExit(
        "filtered_bam_qc.tsv filtered_reads must equal twice filtered_fragments"
    )
if raw_total_reads == 0:
    if filtered_reads:
        raise SystemExit("filtered reads cannot be nonzero when raw total sequences is zero")
    mapq_filtered_fraction = 0.0
else:
    mapq_filtered_fraction = filtered_reads / raw_total_reads
if mapq_filtered_fraction > 1:
    raise SystemExit("MAPQ-filtered reads exceed alignment raw total sequences")

mapped_total = 0
mitochondrial_mapped = 0
with open("idxstats.tsv", encoding="utf-8") as handle:
    for line in handle:
        fields = line.rstrip("\\n").split("\\t")
        if len(fields) < 4 or fields[0] == "*":
            continue
        mapped = int(fields[2])
        mapped_total += mapped
        if fields[0] in {"chrM", "MT", "M"}:
            mitochondrial_mapped += mapped
mitochondrial_percent = (
    100.0 * mitochondrial_mapped / mapped_total if mapped_total else 0.0
)

with open("duplicate_metrics.json", encoding="utf-8") as handle:
    duplicate_json = json.load(handle)
if not isinstance(duplicate_json, dict):
    raise SystemExit("duplicate_metrics.json must contain a JSON object")
duplicate_total = require_non_negative_number(duplicate_json, "DUPLICATE TOTAL")
denominator_key = "EXAMINED" if "EXAMINED" in duplicate_json else "READ"
examined_reads = require_non_negative_number(duplicate_json, denominator_key)
estimated_library_size = (
    require_non_negative_number(duplicate_json, "ESTIMATED LIBRARY SIZE")
    if "ESTIMATED LIBRARY SIZE" in duplicate_json
    else None
)
if duplicate_total > examined_reads:
    raise SystemExit("markdup DUPLICATE TOTAL exceeds examined/read count")
duplicate_percent = (
    100.0 * duplicate_total / examined_reads if examined_reads else 0.0
)

insert_counts = {}
with open("insert_size.tsv", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\\t")
    required_insert_columns = {"record_type", "insert_size", "total_pairs"}
    if reader.fieldnames is None or not required_insert_columns.issubset(reader.fieldnames):
        raise SystemExit(
            "insert_size.tsv must contain record_type, insert_size, and total_pairs"
        )
    for row_number, row in enumerate(reader, start=2):
        if row["record_type"] != "IS":
            raise SystemExit(
                f"insert_size.tsv row {row_number} must have record_type IS"
            )
        try:
            insert_value = int(row["insert_size"])
            pair_count = int(row["total_pairs"])
        except ValueError as error:
            raise SystemExit(
                f"insert_size.tsv row {row_number} values must be integers"
            ) from error
        if insert_value < 0 or pair_count < 0:
            raise SystemExit(
                f"insert_size.tsv row {row_number} values must be non-negative"
            )
        insert_counts[insert_value] = insert_counts.get(insert_value, 0) + pair_count

insert_total_pairs = sum(insert_counts.values())
if insert_total_pairs:
    insert_size_min = min(insert_counts)
    insert_size_max = max(insert_counts)
    insert_size_mean = (
        sum(size * count for size, count in insert_counts.items())
        / insert_total_pairs
    )

    def value_at_rank(rank):
        cumulative = 0
        for size in sorted(insert_counts):
            cumulative += insert_counts[size]
            if rank < cumulative:
                return size
        raise AssertionError("weighted insert-size rank is outside the distribution")

    def weighted_quantile(quantile):
        position = (insert_total_pairs - 1) * quantile
        lower = math.floor(position)
        upper = math.ceil(position)
        lower_value = value_at_rank(lower)
        upper_value = value_at_rank(upper)
        return lower_value + (upper_value - lower_value) * (position - lower)

    insert_size_q25 = weighted_quantile(0.25)
    insert_size_median = weighted_quantile(0.5)
    insert_size_q75 = weighted_quantile(0.75)
else:
    insert_size_min = None
    insert_size_max = None
    insert_size_mean = None
    insert_size_q25 = None
    insert_size_median = None
    insert_size_q75 = None

columns = [
    "sample_id",
    "raw_total_reads",
    "mapped_percent",
    "properly_paired_percent",
    "mapq_filtered_reads",
    "mapq_filtered_fragments",
    "mapq_filtered_fraction",
    "markdup_examined_reads",
    "duplicate_total",
    "duplicate_percent",
    "mitochondrial_percent",
    "estimated_library_size",
    "insert_size_total_pairs",
    "insert_size_min",
    "insert_size_q25",
    "insert_size_mean",
    "insert_size_median",
    "insert_size_q75",
    "insert_size_max",
]
row = {
    "sample_id": sample_id,
    "raw_total_reads": raw_total_reads,
    "mapped_percent": percent_for(" mapped ("),
    "properly_paired_percent": percent_for(" properly paired ("),
    "mapq_filtered_reads": filtered_reads,
    "mapq_filtered_fragments": filtered_fragments,
    "mapq_filtered_fraction": mapq_filtered_fraction,
    "markdup_examined_reads": examined_reads,
    "duplicate_total": duplicate_total,
    "duplicate_percent": duplicate_percent,
    "mitochondrial_percent": mitochondrial_percent,
    "estimated_library_size": estimated_library_size,
    "insert_size_total_pairs": insert_total_pairs,
    "insert_size_min": insert_size_min,
    "insert_size_q25": insert_size_q25,
    "insert_size_mean": insert_size_mean,
    "insert_size_median": insert_size_median,
    "insert_size_q75": insert_size_q75,
    "insert_size_max": insert_size_max,
}
with open(output_name, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\\t")
    writer.writeheader()
    writer.writerow(row)
with open(insert_output_name, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["sample_id", "insert_size", "pair_count"],
        delimiter="\\t",
    )
    writer.writeheader()
    for insert_value in sorted(insert_counts):
        writer.writerow({
            "sample_id": sample_id,
            "insert_size": insert_value,
            "pair_count": insert_counts[insert_value],
        })
PY
    """
}

process MOTIF_QC_CUSTOM {
    tag "${meta.sample_id}"
    label 'process_light'

    conda "${projectDir}/envs/multiqc.yml"
    container 'quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0'

    input:
    tuple val(meta),
        path(motif_qc_tsv, stageAs: 'expected_motif_qc.tsv')

    output:
    tuple val(meta), path("*.motif_qc.tsv"), emit: custom

    script:
    def sampleId = motifQcSampleId(meta)
    def outputName = "${sampleId}.motif_qc.tsv"

    """
    set -euo pipefail
    python - "${sampleId}" "${outputName}" <<'PY'
import csv
import sys

sample_id, output_name = sys.argv[1:]
with open("expected_motif_qc.tsv", encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\\t"))
if not rows or not {"metric", "value"}.issubset(rows[0]):
    raise SystemExit("expected motif QC TSV must contain metric and value columns")
metrics = {}
for row in rows:
    metric = row.get("metric", "")
    if not metric:
        raise SystemExit("expected motif QC TSV contains a blank metric")
    if metric in metrics:
        raise SystemExit(f"expected motif QC TSV contains duplicate metric: {metric}")
    metrics[metric] = row.get("value", "")
status = metrics.get("status", "")
if not status:
    raise SystemExit("expected motif QC TSV is missing status")
row = {
    "sample_id": sample_id,
    "expected_motif_status": status,
    "best_motif_id": metrics.get("best_motif_id", ""),
    "best_adjusted_p_value": metrics.get("best_adjusted_p_value", ""),
}
with open(output_name, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample_id",
            "expected_motif_status",
            "best_motif_id",
            "best_adjusted_p_value",
        ],
        delimiter="\\t",
    )
    writer.writeheader()
    writer.writerow(row)
PY
    """
}

process MULTIQC {
    tag 'run report'
    label 'process_standard'

    conda "${projectDir}/envs/multiqc.yml"
    container 'quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0'

    publishDir "${params.outdir}/reports/multiqc",
        mode: 'copy',
        overwrite: true,
        pattern: 'multiqc_*'
    publishDir "${params.outdir}/reports/summary",
        mode: 'copy',
        overwrite: true,
        pattern: 'combined_target_qc.tsv'
    publishDir "${params.outdir}/reports/summary",
        mode: 'copy',
        overwrite: true,
        pattern: 'peak_enrichment.tsv'
    publishDir "${params.outdir}/reports/summary",
        mode: 'copy',
        overwrite: true,
        pattern: 'matrix_*.tsv'
    publishDir "${params.outdir}/reports/multiqc",
        mode: 'copy',
        overwrite: true,
        pattern: 'observed_vs_null.png'
    publishDir "${params.outdir}/reports/multiqc",
        mode: 'copy',
        overwrite: true,
        pattern: 'matrix_*.png'

    input:
    path fastqc_files, stageAs: 'fastqc??/*'
    path demux_custom_files, stageAs: 'demux??/*'
    path library_custom_files, stageAs: 'library??/*'
    path insert_size_files, stageAs: 'insert??/*'
    path peak_qc_files, stageAs: 'peak_qc??/*'
    path motif_metric_files, stageAs: 'motif??/*'
    path tss_status_files, stageAs: 'tss??/*'
    path enrichment_files, stageAs: 'enrichment??/*'
    val annotation_status

    output:
    path("combined_target_qc.tsv"), emit: combined_summary
    path("peak_enrichment.tsv"), optional: true, emit: enrichment_table
    path("matrix_*.tsv"), optional: true, emit: enrichment_matrices
    path("observed_vs_null.png"), optional: true, emit: enrichment_plot
    path("matrix_*.png"), optional: true, emit: enrichment_heatmaps
    path("multiqc_report.html"), emit: report
    path("multiqc_data"), emit: data
    path("multiqc_custom_content"), emit: custom_content
    path("multiqc_versions.yml"), emit: versions

    script:
    def safeAnnotationStatuses = [
        'skipped_no_annotation',
        'computed_bed',
        'computed_gtf',
    ]
    def annotationStatus = annotation_status.toString()
    if (!(annotationStatus in safeAnnotationStatuses)) {
        throw new IllegalArgumentException(
            "unexpected annotation status: ${annotationStatus}"
        )
    }
    """
    export annotation_status="${annotationStatus}"
    """ +
    '''
    set -euo pipefail
    mkdir -p "multiqc_custom_content"

    python <<'PY'
import csv
import glob
import os
import shutil
from pathlib import Path

custom_dir = Path("multiqc_custom_content")


def read_table(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def combine_tables(pattern, identity):
    combined = []
    seen = set()
    for path in sorted(glob.glob(pattern)):
        for row in read_table(path):
            key = row.get(identity, "")
            if not key:
                raise SystemExit(f"{path}: missing {identity}")
            if key in seen:
                raise SystemExit(f"duplicate {identity} in custom QC inputs: {key}")
            seen.add(key)
            combined.append(row)
    return sorted(combined, key=lambda row: row[identity])


def append_yaml_sections(handle, sections):
    handle.write("custom_data:\n")
    for key, metadata in sections:
        handle.write(f"  {key}:\n")
        for field, value in metadata.items():
            escaped = value.replace('"', '\\"')
            handle.write(f'    {field}: "{escaped}"\n')
    handle.write("sp:\n")
    for key, metadata in sections:
        handle.write(f"  {key}:\n")
        handle.write(f'    fn: "{metadata["filename"]}"\n')
    handle.write("ignore_images: false\n")


demux_columns = [
    "library_id",
    "total_read_pairs",
    "assigned_read_pairs",
    "ambiguous_read_pairs",
    "unassigned_read_pairs",
    "assigned_fraction",
    "ambiguous_fraction",
    "unassigned_fraction",
]
demux_rows = combine_tables("demux*/*.demultiplex_qc.tsv", "library_id")
write_table(
    custom_dir / "nanocut_demultiplex_mqc.tsv",
    demux_columns,
    demux_rows,
)

library_columns = [
    "sample_id",
    "raw_total_reads",
    "mapped_percent",
    "properly_paired_percent",
    "mapq_filtered_reads",
    "mapq_filtered_fragments",
    "mapq_filtered_fraction",
    "markdup_examined_reads",
    "duplicate_total",
    "duplicate_percent",
    "mitochondrial_percent",
    "estimated_library_size",
    "insert_size_total_pairs",
    "insert_size_min",
    "insert_size_q25",
    "insert_size_mean",
    "insert_size_median",
    "insert_size_q75",
    "insert_size_max",
]
library_rows = combine_tables("library*/*.library_qc.tsv", "sample_id")
write_table(
    custom_dir / "nanocut_library_qc_mqc.tsv",
    library_columns,
    library_rows,
)

insert_by_sample = {}
for path in sorted(glob.glob("insert*/*.insert_size_distribution.tsv")):
    rows = read_table(path)
    sample_ids = {row.get("sample_id", "") for row in rows}
    if not rows:
        sample_id = Path(path).name.removesuffix(".insert_size_distribution.tsv")
        sample_ids = {sample_id}
    if len(sample_ids) != 1 or "" in sample_ids:
        raise SystemExit(f"{path}: insert-size rows must contain one sample_id")
    sample_id = next(iter(sample_ids))
    if sample_id in insert_by_sample:
        raise SystemExit(f"duplicate insert-size sample_id: {sample_id}")
    counts = {}
    for row in rows:
        try:
            insert_size = int(row["insert_size"])
            pair_count = int(row["pair_count"])
        except (KeyError, ValueError) as error:
            raise SystemExit(
                f"{path}: insert_size and pair_count must be integers"
            ) from error
        if insert_size < 0 or pair_count < 0:
            raise SystemExit(
                f"{path}: insert_size and pair_count must be non-negative"
            )
        if insert_size in counts:
            raise SystemExit(
                f"{path}: duplicate insert_size {insert_size} for {sample_id}"
            )
        counts[insert_size] = pair_count
    insert_by_sample[sample_id] = counts
insert_samples = sorted(insert_by_sample)
insert_sizes = sorted({
    insert_size
    for counts in insert_by_sample.values()
    for insert_size in counts
})
insert_rows = [
    {
        "insert_size": insert_size,
        **{
            sample_id: insert_by_sample[sample_id].get(insert_size, 0)
            for sample_id in insert_samples
        },
    }
    for insert_size in insert_sizes
]
write_table(
    custom_dir / "nanocut_insert_size_mqc.tsv",
    ["insert_size", *insert_samples],
    insert_rows,
)

peak_rows = []
seen_peak_samples = set()
for path in sorted(glob.glob("peak_qc*/*.peak_qc.tsv")):
    metrics = {}
    for row in read_table(path):
        metric = row.get("metric")
        if metric:
            metrics[metric] = row.get("value", "")
    sample_id = metrics.get("sample_id", "")
    if not sample_id:
        raise SystemExit(f"{path}: peak-QC TSV has no sample_id metric")
    if sample_id in seen_peak_samples:
        raise SystemExit(f"duplicate peak-QC sample_id: {sample_id}")
    seen_peak_samples.add(sample_id)
    peak_rows.append(metrics)
peak_rows.sort(key=lambda row: row["sample_id"])

priority = [
    "sample_id",
    "peak_count",
    "total_covered_bases",
    "total_fragments",
    "fragments_in_peaks",
    "frip",
]
all_peak_metrics = set().union(*(row.keys() for row in peak_rows)) if peak_rows else set()
combined_columns = priority + sorted(all_peak_metrics.difference(priority))
write_table("combined_target_qc.tsv", combined_columns, peak_rows)
write_table(
    custom_dir / "nanocut_peak_qc_mqc.tsv",
    ["sample_id", "frip", "peak_count", "total_covered_bases"],
    peak_rows,
)

motif_by_sample = {}
for path in sorted(glob.glob("motif*/*.motif_qc.tsv")):
    for row in read_table(path):
        sample_id = row.get("sample_id", "")
        if not sample_id:
            raise SystemExit(f"{path}: motif custom TSV is missing sample_id")
        if sample_id in motif_by_sample:
            raise SystemExit(f"duplicate motif sample_id: {sample_id}")
        motif_by_sample[sample_id] = row
motif_rows = []
for peak_row in peak_rows:
    sample_id = peak_row["sample_id"]
    motif_rows.append(
        motif_by_sample.get(
            sample_id,
            {
                "sample_id": sample_id,
                "expected_motif_status": "not_run",
                "best_motif_id": "",
                "best_adjusted_p_value": "",
            },
        )
    )
write_table(
    custom_dir / "nanocut_motif_qc_mqc.tsv",
    [
        "sample_id",
        "expected_motif_status",
        "best_motif_id",
        "best_adjusted_p_value",
    ],
    motif_rows,
)

tss_rows = []
for path in sorted(glob.glob("tss*/*.tss_status.tsv")):
    tss_rows.extend(read_table(path))
if not tss_rows:
    tss_rows = [{
        "sample_id": "<run>",
        "annotation_mode": "none",
        "status": os.environ.get("annotation_status", "skipped_no_annotation"),
    }]
write_table(
    custom_dir / "nanocut_tss_qc_mqc.tsv",
    ["sample_id", "annotation_mode", "status"],
    tss_rows,
)

enrichment_columns = [
    "comparison_id",
    "foreground_id",
    "foreground_tf",
    "reference_id",
    "reference_type",
    "reference_tf",
    "background_model",
    "foreground_peak_count",
    "reference_peak_count",
    "observed_overlap_count",
    "null_mean_overlap",
    "null_sd_overlap",
    "enrichment_ratio",
    "empirical_p_value",
    "permutations_requested",
    "permutations_succeeded",
    "seed",
    "status",
]
enrichment_rows = []
seen_enrichment_comparisons = set()
copied_enrichment = set()
heatmap_models = [
    (
        "random",
        "Fully random genomic intervals.",
        "matrix_random.png",
        "matrix_random.tsv",
        "nanocut_peak_enrichment_random_heatmap_mqc.png",
    ),
    (
        "length_matched",
        "Random intervals matched only on peak length.",
        "matrix_length_matched.png",
        "matrix_length_matched.tsv",
        "nanocut_peak_enrichment_length_matched_heatmap_mqc.png",
    ),
    (
        "gc_matched",
        "Random intervals matched only on GC content.",
        "matrix_gc_matched.png",
        "matrix_gc_matched.tsv",
        "nanocut_peak_enrichment_gc_matched_heatmap_mqc.png",
    ),
    (
        "length_gc_matched",
        "Random intervals matched on both length and GC content.",
        "matrix_length_gc_matched.png",
        "matrix_length_gc_matched.tsv",
        "nanocut_peak_enrichment_length_gc_matched_heatmap_mqc.png",
    ),
]


def register_enrichment_file(path):
    destination = Path(path.name)
    if destination.name in copied_enrichment:
        return
    shutil.copyfile(path, destination)
    copied_enrichment.add(destination.name)
    if destination.name == "peak_enrichment.tsv":
        for row in read_table(path):
            if (
                row.get("foreground_id")
                and row.get("reference_id")
                and row.get("background_model")
                and row.get("status") == "ok"
            ):
                comparison_id = (
                    f"{row['foreground_id']}|{row['reference_id']}|{row['background_model']}"
                )
                if comparison_id in seen_enrichment_comparisons:
                    raise SystemExit(
                        f"duplicate enrichment comparison_id: {comparison_id}"
                    )
                seen_enrichment_comparisons.add(comparison_id)
                enrichment_rows.append({"comparison_id": comparison_id, **row})


for path_text in sorted(glob.glob("enrichment*/*")):
    path = Path(path_text)
    if path.is_dir():
        for nested in sorted(path.iterdir()):
            if nested.is_file():
                register_enrichment_file(nested)
    elif path.is_file():
        register_enrichment_file(path)

if enrichment_rows:
    for model, description, heatmap_png, matrix_tsv, heatmap_asset in heatmap_models:
        heatmap_path = Path(heatmap_png)
        matrix_path = Path(matrix_tsv)
        if not heatmap_path.is_file() or not matrix_path.is_file():
            raise SystemExit(
                f"missing enrichment heatmap assets for {model}: "
                f"{heatmap_png}, {matrix_tsv}"
            )
        shutil.copyfile(heatmap_path, custom_dir / heatmap_asset)
    write_table(
        custom_dir / "nanocut_peak_enrichment_mqc.tsv",
        enrichment_columns,
        enrichment_rows,
    )
    with open(
        custom_dir / "nanocut_peak_enrichment_overview_mqc.md",
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write('id: "nanocut_peak_enrichment_overview"\n')
        handle.write('section_name: "Nano-CUT&Tag peak enrichment"\n')
        handle.write('description: "Peak-overlap enrichment against four null models."\n')
        handle.write("---\n")
        handle.write(
            "Peak enrichment is reported in "
            "[`peak_enrichment.tsv`](../summary/peak_enrichment.tsv) "
            "with one row per foreground/reference/null-model comparison.\n\n"
        )
        handle.write(
            "Null models and heatmaps:\n"
        )
        for model, description, heatmap_png, matrix_tsv, heatmap_asset in heatmap_models:
            handle.write(
                f"- `{model}` — {description} "
                f"[heatmap]({heatmap_png}), "
                f"[matrix](../summary/{matrix_tsv})\n"
            )
        handle.write(
            "\nThe observed/null overview plot is available as "
            "[`observed_vs_null.png`](observed_vs_null.png).\n"
        )

config_sections = [
    (
        "nanocut_demultiplex",
        {
            "section_name": "Nano-CUT&Tag demultiplexing",
            "description": "Physical-library barcode assignment metrics.",
            "plot_type": "table",
            "file_format": "tsv",
            "filename": "nanocut_demultiplex_mqc.tsv",
        },
    ),
    (
        "nanocut_library_qc",
        {
            "section_name": "Nano-CUT&Tag library QC",
            "description": "Alignment, pairing, duplicate, mitochondrial, and complexity metrics.",
            "plot_type": "table",
            "file_format": "tsv",
            "filename": "nanocut_library_qc_mqc.tsv",
        },
    ),
    (
        "nanocut_insert_size",
        {
            "section_name": "Nano-CUT&Tag insert-size distribution",
            "description": "Paired-fragment insert-size counts from SAMtools stats.",
            "plot_type": "linegraph",
            "file_format": "tsv",
            "filename": "nanocut_insert_size_mqc.tsv",
        },
    ),
    (
        "nanocut_peak_qc",
        {
            "section_name": "Nano-CUT&Tag broad-peak QC",
            "description": "Final broad-peak counts and fragment-based FRiP.",
            "plot_type": "table",
            "file_format": "tsv",
            "filename": "nanocut_peak_qc_mqc.tsv",
        },
    ),
    (
        "nanocut_motif_qc",
        {
            "section_name": "Nano-CUT&Tag expected motif QC",
            "description": "Expected-motif status; not_run is retained until motif results are supplied.",
            "plot_type": "table",
            "file_format": "tsv",
            "filename": "nanocut_motif_qc_mqc.tsv",
        },
    ),
    (
        "nanocut_tss_qc",
        {
            "section_name": "Nano-CUT&Tag TSS enrichment",
            "description": "Optional strand-aware TSS enrichment status.",
            "plot_type": "table",
            "file_format": "tsv",
            "filename": "nanocut_tss_qc_mqc.tsv",
        },
    ),
]
if enrichment_rows:
    config_sections.append(
        (
            "nanocut_peak_enrichment",
            {
                "section_name": "Nano-CUT&Tag peak enrichment",
                "description": "Foreground/reference overlap enrichment across the random, length_matched, gc_matched, and length_gc_matched null models.",
                "plot_type": "table",
                "file_format": "tsv",
                "filename": "nanocut_peak_enrichment_mqc.tsv",
            },
        )
    )
    for model, description, heatmap_png, matrix_tsv, heatmap_asset in heatmap_models:
        config_sections.append(
            (
                f"nanocut_peak_enrichment_{model}_heatmap",
                {
                    "section_name": f"Peak enrichment heatmap: {model}",
                    "description": description,
                    "plot_type": "image",
                    "file_format": "png",
                    "filename": heatmap_asset,
                },
            )
        )

with open("multiqc_config.yml", "w", encoding="utf-8") as handle:
    append_yaml_sections(handle, config_sections)
PY

    cp "multiqc_config.yml" "multiqc_custom_content/multiqc_config.yml"

    multiqc \
        --force \
        --config "multiqc_config.yml" \
        --outdir . \
        .

    printf 'MULTIQC:\n  multiqc: ' > "multiqc_versions.yml"
    multiqc --version 2>&1 | sed -n '1p' >> "multiqc_versions.yml"
    '''
}
