process BAMCOVERAGE {
    tag "${meta.sample_id}"
    label 'process_heavy'

    conda "${projectDir}/envs/deeptools.yml"
    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'

    publishDir "${params.outdir}/coverage/${meta.sample_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(filtered_bam, stageAs: 'coverage.filtered.bam'),
        path(filtered_bai, stageAs: 'coverage.filtered.bam.bai')
    val min_mapq

    output:
    tuple val(meta), path("coverage.RPKM.bw"), emit: bigwig
    tuple val(meta), path("bamcoverage_versions.yml"), emit: versions

    script:
    def mapqText = min_mapq.toString()
    if (!(mapqText ==~ /[0-9]+/)) {
        throw new IllegalArgumentException(
            "min_mapq must be a non-negative integer, got ${min_mapq}"
        )
    }
    def mapq
    try {
        mapq = Integer.parseInt(mapqText)
    } catch (NumberFormatException ignored) {
        throw new IllegalArgumentException(
            "min_mapq is outside the supported integer range: ${min_mapq}"
        )
    }

    """
    set -euo pipefail
    bamCoverage \
        --bam "coverage.filtered.bam" \
        --outFileName "coverage.RPKM.bw" \
        --numberOfProcessors "${task.cpus}" \
        --minMappingQuality "${mapq}" \
        --binSize 50 \
        --centerReads \
        --smoothLength 250 \
        --normalizeUsing RPKM \
        --ignoreDuplicates \
        --extendReads

    printf 'BAMCOVERAGE:\\n  deeptools: ' > "bamcoverage_versions.yml"
    bamCoverage --version 2>&1 | sed -n '1p' >> "bamcoverage_versions.yml"
    """
}
