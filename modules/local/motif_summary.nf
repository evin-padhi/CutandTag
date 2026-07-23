def motifSummarySampleId(meta) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    def sampleId = meta instanceof Map ? meta.sample_id?.toString() : null
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "motif-summary sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    sampleId
}

def motifSummaryExpectedPattern(meta) {
    def expectedMotif = meta instanceof Map
        ? meta.expected_motif?.toString()
        : null
    if (expectedMotif == null || expectedMotif.trim().isEmpty()) {
        throw new IllegalArgumentException(
            "motif-summary expected_motif must not be blank"
        )
    }
    expectedMotif
}

def motifSummaryWindow(rawWindow) {
    def text = rawWindow?.toString()
    if (text == null || !(text ==~ /[0-9]+/)) {
        throw new IllegalArgumentException(
            "motif-summary window must be a positive integer, got ${text}"
        )
    }
    def value = new BigInteger(text)
    if (value.signum() <= 0 || value > Integer.MAX_VALUE) {
        throw new IllegalArgumentException(
            "motif-summary window must be a positive integer, got ${text}"
        )
    }
    value.toString()
}

process MOTIF_SUMMARY {
    tag "${meta.sample_id}"
    label 'process_light'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/motifs/${meta.sample_id}/expected_motif_qc",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(ame_results, stageAs: 'ame'),
        path(fimo_results, stageAs: 'fimo')
    val motif_window

    output:
    tuple val(meta),
        path("expected_motif_qc.json"),
        path("expected_motif_qc.tsv"),
        path("expected_motif_qc.motif_hit_positions.tsv"),
        emit: qc
    tuple val(meta), path("motif_summary_versions.yml"), emit: versions

    script:
    motifSummarySampleId(meta)
    def expectedMotif = motifSummaryExpectedPattern(meta)
    def encodedExpectedMotif = java.util.Base64.getEncoder().encodeToString(
        expectedMotif.getBytes('UTF-8')
    )
    def window = motifSummaryWindow(motif_window)

    """
    set -euo pipefail
    printf '%s' "${encodedExpectedMotif}" > "expected_motif.b64"
    expected_motif=\$(base64 --decode "expected_motif.b64")

    motif_qc.py \
        --ame "ame/ame.tsv" \
        --fimo "fimo/fimo.tsv" \
        --expected-motif "\${expected_motif}" \
        --window "${window}" \
        --output-prefix "expected_motif_qc"

    printf 'MOTIF_SUMMARY:\\n  python: ' > "motif_summary_versions.yml"
    python --version 2>&1 >> "motif_summary_versions.yml"
    printf '  motif_qc.py: repository\\n' >> "motif_summary_versions.yml"
    """
}
