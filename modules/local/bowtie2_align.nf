process BOWTIE2_ALIGN {
    tag "${meta.sample_id}"
    label 'process_heavy'

    conda "${projectDir}/envs/bowtie2.yml"
    container 'quay.io/biocontainers/bowtie2:2.5.4--he20e202_2'

    publishDir "${params.outdir}/alignment/${meta.sample_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(r1, stageAs: 'reads_R1.fastq.gz'),
        path(r2, stageAs: 'reads_R2.fastq.gz'),
        val(reference_meta),
        path(index_files, stageAs: 'bt2_index/*')

    output:
    tuple val(meta), path("alignment.sam"), emit: sam
    tuple val(meta), path("bowtie2.summary.txt"), emit: summary
    tuple val(meta), path("bowtie2_align_versions.yml"), emit: versions

    script:
    """
    set -euo pipefail
    bowtie2 \
        --threads "${task.cpus}" \
        -x "bt2_index/reference" \
        -1 "reads_R1.fastq.gz" \
        -2 "reads_R2.fastq.gz" \
        -S "alignment.sam" \
        2> "bowtie2.summary.txt"

    printf 'BOWTIE2_ALIGN:\\n  bowtie2: ' > "bowtie2_align_versions.yml"
    bowtie2 --version 2>&1 | sed -n '1p' >> "bowtie2_align_versions.yml"
    """
}
