process BIGWIG_SUBTRACT_IGG {
    tag "${meta.sample_id} minus ${control_meta.sample_id}"
    label 'process_heavy'

    conda "${projectDir}/envs/deeptools.yml"
    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'

    publishDir path: {
        "${params.outdir}/coverage/igg_subtracted/${meta.sample_id}"
    },
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta), val(control_meta),
        path(target_bigwig, stageAs: 'target.RPKM.bw'),
        path(control_bigwig, stageAs: 'matched_igg.RPKM.bw')

    output:
    tuple val(meta), val(control_meta),
        path("${meta.sample_id}.RPKM_minus_IgG.bw"),
        emit: bigwig
    tuple val(meta), path("bigwig_subtract_igg_versions.yml"), emit: versions

    script:
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    def sampleId = meta.sample_id?.toString()
    def controlId = control_meta.sample_id?.toString()
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "target sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    if (controlId == null || !controlId.matches(safeId)) {
        throw new IllegalArgumentException(
            "control sample_id must match ${safeId}, got ${controlId}"
        )
    }
    if (meta.control_id?.toString() != controlId) {
        throw new IllegalArgumentException(
            "target ${sampleId} control_id does not match ${controlId}"
        )
    }
    if (meta.is_control || !control_meta.is_control ||
        control_meta.assay_target?.toString() != 'IgG') {
        throw new IllegalArgumentException(
            "${sampleId} must pair with a matched IgG control"
        )
    }
    def outputName = "${sampleId}.RPKM_minus_IgG.bw"

    """
    set -euo pipefail
    printf 'Subtracting matched IgG RPKM signal: %s - %s\\n' \\
        "${sampleId}" "${controlId}"
    bigwigCompare \\
        --bigwig1 "target.RPKM.bw" \\
        --bigwig2 "matched_igg.RPKM.bw" \\
        --operation subtract \\
        --binSize 50 \\
        --skipZeroOverZero \\
        --numberOfProcessors "${task.cpus}" \\
        --outFileName "${outputName}"

    printf 'BIGWIG_SUBTRACT_IGG:\\n  deeptools: ' \\
        > "bigwig_subtract_igg_versions.yml"
    bigwigCompare --version 2>&1 | sed -n '1p' \\
        >> "bigwig_subtract_igg_versions.yml"
    """
}
