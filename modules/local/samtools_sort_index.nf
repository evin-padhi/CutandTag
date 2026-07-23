process SAMTOOLS_SORT_INDEX {
    tag "${meta.sample_id}"
    label 'process_heavy'

    conda "${projectDir}/envs/samtools.yml"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    publishDir "${params.outdir}/alignment/${meta.sample_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta), path(sam, stageAs: 'source.alignment.sam')

    output:
    tuple val(meta), path("analysis.bam"), path("analysis.bam.bai"), emit: bam
    tuple val(meta), path("samtools_sort_index_versions.yml"), emit: versions

    script:
    """
    set -euo pipefail
    samtools sort \
        -@ "${task.cpus}" \
        -o "analysis.bam" \
        "source.alignment.sam"
    samtools quickcheck -v "analysis.bam"
    samtools index \
        -@ "${task.cpus}" \
        "analysis.bam" \
        "analysis.bam.bai"

    printf 'SAMTOOLS_SORT_INDEX:\\n  samtools: ' > "samtools_sort_index_versions.yml"
    samtools --version 2>&1 | sed -n '1p' >> "samtools_sort_index_versions.yml"
    """
}
