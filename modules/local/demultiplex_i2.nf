process DEMULTIPLEX_I2 {
    tag "${meta.library_id}"
    label 'process_standard'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/demultiplex/${meta.library_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta), path(r1, stageAs: 'input_R1.fastq.gz'), path(r2, stageAs: 'input_R2.fastq.gz'), path(i2, stageAs: 'input_I2.fastq.gz'), val(samples)
    val max_mismatches
    val allow_empty

    output:
    tuple val(meta), val(samples), path("demux/*_R1.fastq.gz"), path("demux/*_R2.fastq.gz"), emit: reads
    tuple val(meta), path("demultiplex.metrics.json"), path("demultiplex.metrics.tsv"), emit: metrics
    path "demultiplex_i2_versions.yml", emit: versions

    script:
    def mismatchText = max_mismatches.toString()
    if (!(mismatchText ==~ /[0-9]+/)) {
        throw new IllegalArgumentException(
            "max_mismatches must be a non-negative integer, got ${max_mismatches}"
        )
    }
    def mismatchCount
    try {
        mismatchCount = Integer.parseInt(mismatchText)
    } catch (NumberFormatException ignored) {
        throw new IllegalArgumentException(
            "max_mismatches is outside the supported integer range: ${max_mismatches}"
        )
    }
    def sampleArgs = samples.collect { sample ->
        def token = "${sample.sample_id}=${sample.barcode}"
        def escaped = token.replace("'", "'\"'\"'")
        "--sample '${escaped}'"
    }.join(' ')
    def allowEmptyFlag = allow_empty.toString().toBoolean() ? '--allow-empty' : ''

    """
    mkdir -p demux
    demultiplex_i2.py \
        --r1 "input_R1.fastq.gz" \
        --r2 "input_R2.fastq.gz" \
        --i2 "input_I2.fastq.gz" \
        ${sampleArgs} \
        --max-mismatches "${mismatchCount}" \
        --outdir "demux" \
        --metrics "demultiplex.metrics.json" \
        ${allowEmptyFlag}

    printf 'DEMULTIPLEX_I2:\\n  python: ' > demultiplex_i2_versions.yml
    python --version 2>&1 >> demultiplex_i2_versions.yml
    printf '  demultiplex_i2.py: repository\\n' >> demultiplex_i2_versions.yml
    """
}
