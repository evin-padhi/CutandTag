process SAMTOOLS_FILTER {
    tag "${meta.sample_id}"
    label 'process_standard'

    conda "${projectDir}/envs/samtools.yml"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    publishDir "${params.outdir}/alignment/${meta.sample_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(analysis_bam, stageAs: 'source.analysis.bam'),
        path(analysis_bai, stageAs: 'source.analysis.bam.bai')
    val min_mapq

    output:
    tuple val(meta), path("filtered.bam"), path("filtered.bam.bai"), emit: bam
    tuple val(meta), path("samtools_filter_versions.yml"), emit: versions

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
    samtools view \
        -@ "${task.cpus}" \
        -b \
        -f 2 \
        -F 2304 \
        -q "${mapq}" \
        -o "filtered.bam" \
        "source.analysis.bam"
    samtools quickcheck -v "filtered.bam"
    samtools index \
        -@ "${task.cpus}" \
        "filtered.bam" \
        "filtered.bam.bai"

    printf 'SAMTOOLS_FILTER:\\n  samtools: ' > "samtools_filter_versions.yml"
    samtools --version 2>&1 | sed -n '1p' >> "samtools_filter_versions.yml"
    """
}
