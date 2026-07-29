process FASTQC {
    tag "${meta.sample_id}"
    label 'process_standard'

    conda "${projectDir}/envs/fastqc.yml"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    publishDir "${params.outdir}/fastqc/${meta.sample_id}",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename.tokenize('/').last() }

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("fastqc/*_fastqc.html"), path("fastqc/*_fastqc.zip"), emit: reports
    path "fastqc_versions.yml", emit: versions

    script:
    """
    mkdir -p fastqc
    fastqc \
        --threads ${task.cpus} \
        --outdir fastqc \
        "${r1}" \
        "${r2}"

    printf 'FASTQC:\\n  fastqc: ' > fastqc_versions.yml
    fastqc --version 2>&1 | sed 's/^FastQC v//' >> fastqc_versions.yml
    """
}
