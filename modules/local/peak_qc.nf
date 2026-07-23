def peakQcSampleId(meta) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    if (!(meta instanceof Map)) {
        throw new IllegalArgumentException(
            "peak-QC metadata must be a map"
        )
    }
    def sampleId = meta.sample_id?.toString()
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "peak-QC sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    sampleId
}

process PEAK_QC {
    tag "${meta.sample_id}"
    label 'process_light'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/qc/peaks/${meta.sample_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(fragments, stageAs: 'fragments.bedpe'),
        path(final_broad_peaks, stageAs: 'final.broadPeak')

    output:
    tuple val(meta),
        path("*.peak_qc.json"),
        path("*.peak_qc.tsv"),
        path("*.peak_qc.width_histogram.tsv"),
        path("*.peak_qc.fragments_per_peak.tsv"),
        emit: qc
    tuple val(meta), path("peak_qc_versions.yml"), emit: versions

    script:
    def sampleId = peakQcSampleId(meta)
    def outputStem = "${sampleId}.peak_qc"

    """
    set -euo pipefail
    peak_qc.py \
        --fragments "fragments.bedpe" \
        --peaks "final.broadPeak" \
        --sample-id "${sampleId}" \
        --output-prefix "${outputStem}"

    printf 'PEAK_QC:\\n  python: ' > "peak_qc_versions.yml"
    python --version 2>&1 >> "peak_qc_versions.yml"
    printf '  peak_qc.py: repository\\n' >> "peak_qc_versions.yml"
    """
}
