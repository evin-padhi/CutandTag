process QC_DASHBOARD {
    tag 'consolidated QC dashboard'
    label 'process_standard'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/reports/qc_dashboard",
        mode: 'copy',
        overwrite: true,
        pattern: '{qc_dashboard.html,qc_summary.tsv,qc_summary.json,top_motifs.tsv,tss_profiles.tsv}'

    input:
    val encoded_metadata
    path demultiplex_jsons,
        stageAs: 'dashboard_inputs/demultiplex/demux??/*'
    path library_metrics,
        stageAs: 'dashboard_inputs/library/library??/*'
    path insert_sizes,
        stageAs: 'dashboard_inputs/insert/insert??/*'
    path peak_metrics,
        stageAs: 'dashboard_inputs/peak/metrics??/*'
    path peak_width_histograms,
        stageAs: 'dashboard_inputs/peak/widths??/*'
    path tss_profile_data,
        stageAs: 'dashboard_inputs/tss/profiles??/*'
    path tss_statuses,
        stageAs: 'dashboard_inputs/tss/statuses??/*'
    path expected_motif_metrics,
        stageAs: 'dashboard_inputs/motif/metrics??/*'
    path ame_results,
        stageAs: 'dashboard_inputs/ame/results??/*'
    val ame_result_sample_ids
    path ame_statuses,
        stageAs: 'dashboard_inputs/ame_status/statuses??/*'
    val ame_status_sample_ids
    val annotation_status

    output:
    path "qc_dashboard.html", emit: report
    path "qc_summary.tsv", emit: summary_tsv
    path "qc_summary.json", emit: summary_json
    path "top_motifs.tsv", emit: top_motifs
    path "tss_profiles.tsv", emit: tss_profiles
    path "qc_dashboard_versions.yml", emit: versions

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
    def encodedMetadata = encoded_metadata.toString()
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    def resultSampleIds = ame_result_sample_ids.collect { it.toString() }
    def statusSampleIds = ame_status_sample_ids.collect { it.toString() }
    if (resultSampleIds.any { !it.matches(safeId) }) {
        throw new IllegalArgumentException(
            "AME result sample_ids contain an unsafe value"
        )
    }
    if (statusSampleIds.any { !it.matches(safeId) }) {
        throw new IllegalArgumentException(
            "AME status sample_ids contain an unsafe value"
        )
    }
    def resultSampleIdsJson = groovy.json.JsonOutput.toJson(resultSampleIds)
    def statusSampleIdsJson = groovy.json.JsonOutput.toJson(statusSampleIds)

    """
    set -euo pipefail
    mkdir -p \
        "dashboard_inputs/demultiplex" \
        "dashboard_inputs/library" \
        "dashboard_inputs/insert" \
        "dashboard_inputs/peak" \
        "dashboard_inputs/tss" \
        "dashboard_inputs/motif" \
        "dashboard_inputs/ame" \
        "dashboard_inputs/ame_status" \
        "dashboard_inputs/normalized_motif" \
        "dashboard_inputs/normalized_ame"

    printf '%s' '${encodedMetadata}' | python -c \
        'import base64, sys; sys.stdout.buffer.write(base64.b64decode(sys.stdin.buffer.read(), validate=True))' \
        > "sample_metadata.json"

    python <<'PY'
import csv
import json
import shutil
from pathlib import Path

with open("sample_metadata.json", encoding="utf-8") as handle:
    metadata = json.load(handle)
target_ids = sorted(
    record["sample_id"]
    for record in metadata
    if not record["is_control"]
)
result_sample_ids = json.loads(r'''${resultSampleIdsJson}''')
status_sample_ids = json.loads(r'''${statusSampleIdsJson}''')
if len(result_sample_ids) != len(set(result_sample_ids)):
    raise SystemExit("duplicate AME result sample_id")
if len(status_sample_ids) != len(set(status_sample_ids)):
    raise SystemExit("duplicate AME status sample_id")
if bool(result_sample_ids) != bool(status_sample_ids):
    raise SystemExit("AME result and status identities must both be empty or non-empty")
if result_sample_ids and sorted(result_sample_ids) != target_ids:
    raise SystemExit("AME result sample_ids do not match target metadata")
if status_sample_ids and sorted(status_sample_ids) != target_ids:
    raise SystemExit("AME status sample_ids do not match target metadata")

ame_root = Path("dashboard_inputs/ame")
ame_tables = []
for staged in sorted(ame_root.glob("results*/*")):
    candidate = staged / "ame.tsv" if staged.is_dir() else staged
    if candidate.name == "ame.tsv" and candidate.is_file():
        ame_tables.append(candidate)

status_tables = sorted(
    path
    for path in Path("dashboard_inputs/ame_status").glob("statuses*/*")
    if path.name == "ame_status.tsv" and path.is_file()
)
if len(ame_tables) != len(result_sample_ids):
    raise SystemExit(
        f"AME result count {len(ame_tables)} does not match "
        f"sample_id count {len(result_sample_ids)}"
    )
if len(status_tables) != len(status_sample_ids):
    raise SystemExit(
        f"AME status count {len(status_tables)} does not match "
        f"sample_id count {len(status_sample_ids)}"
    )

normalized_ame = Path("dashboard_inputs/normalized_ame")
for sample_id, source in zip(result_sample_ids, ame_tables):
    destination = normalized_ame / sample_id / "ame"
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination / "ame.tsv")

ame_status_by_sample = {}
for sample_id, path in zip(status_sample_ids, status_tables):
    status_bytes = path.read_bytes()
    allowed_statuses = {
        b"status\\tcomputed\\n": "computed",
        b"status\\tno_peaks\\n": "no_peaks",
    }
    if status_bytes not in allowed_statuses:
        raise SystemExit(
            f"{path}: AME status must be exactly status followed by "
            "computed or no_peaks"
        )
    ame_status_by_sample[sample_id] = allowed_statuses[status_bytes]

normalized_motif = Path("dashboard_inputs/normalized_motif")
for source in sorted(Path("dashboard_inputs/motif").glob("metrics*/*")):
    if not source.name.endswith(".motif_qc.tsv") or not source.is_file():
        continue
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        columns = list(reader.fieldnames or [])
    if len(rows) != 1 or not rows[0].get("sample_id", "").strip():
        raise SystemExit(f"{source}: motif metrics must contain one sample_id row")
    sample_id = rows[0]["sample_id"].strip()
    if sample_id in ame_status_by_sample:
        rows[0]["ame_status"] = ame_status_by_sample[sample_id]
        if "ame_status" not in columns:
            columns.append("ame_status")
    destination = normalized_motif / source.name
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
PY

    qc_dashboard.py \
        --metadata "sample_metadata.json" \
        --demux-dir "dashboard_inputs/demultiplex" \
        --library-dir "dashboard_inputs/library" \
        --insert-dir "dashboard_inputs/insert" \
        --peak-dir "dashboard_inputs/peak" \
        --tss-dir "dashboard_inputs/tss" \
        --motif-dir "dashboard_inputs/normalized_motif" \
        --ame-dir "dashboard_inputs/normalized_ame" \
        --annotation-status "${annotationStatus}" \
        --outdir "."

    cat > "qc_dashboard_versions.yml" <<'EOF'
QC_DASHBOARD:
  python: Python 3.12.3
  qc_dashboard.py: repository
EOF
    """
}
