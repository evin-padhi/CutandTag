process DIRECT_INPUT_METRICS {
    tag "${meta.library_id}"
    label 'process_light'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir path: { "${params.outdir}/demultiplex/${meta.library_id}" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename.tokenize('/').last() }

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path('demultiplex.metrics.json'),
        path('demultiplex.metrics.tsv'), emit: metrics
    path 'direct_input_versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    printf 'DIRECT_INPUT_METRICS: creating status record for already split sample %s\\n' '${meta.sample_id}'
    python - '${meta.library_id}' '${meta.sample_id}' <<'PY'
import csv
import json
import sys

library_id, sample_id = sys.argv[1:]
payload = {
    "processing_mode": "already_split",
    "library_id": library_id,
    "sample_id": sample_id,
    "total_reads": None,
    "assigned_reads": None,
    "ambiguous_reads": None,
    "unassigned_reads": None,
    "assigned_fraction": None,
    "ambiguous_fraction": None,
    "unassigned_fraction": None,
}
with open("demultiplex.metrics.json", "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2)
    handle.write("\\n")
with open("demultiplex.metrics.tsv", "w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t")
    writer.writerow(["metric", "value"])
    writer.writerow(["processing_mode", "already_split"])
    writer.writerow(["sample_id", sample_id])
PY

    printf 'DIRECT_INPUT_METRICS:\\n  python: ' > direct_input_versions.yml
    python --version 2>&1 >> direct_input_versions.yml
    printf '  direct input metrics: repository\\n' >> direct_input_versions.yml
    """
}
