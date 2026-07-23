def stremeSampleId(meta) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    def sampleId = meta instanceof Map ? meta.sample_id?.toString() : null
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "STREME sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    sampleId
}

process STREME {
    tag "${meta.sample_id}"
    label 'process_heavy'

    conda "${projectDir}/envs/meme.yml"
    container 'quay.io/biocontainers/meme:5.5.7--pl5321h1ca524f_3'

    publishDir "${params.outdir}/motifs/${meta.sample_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(foreground, stageAs: 'foreground.fa'),
        path(background, stageAs: 'background.fa'),
        path(sequence_status, stageAs: 'motif_sequence_status.tsv')

    output:
    tuple val(meta), path("streme"), emit: results
    tuple val(meta), path("streme_status.tsv"), emit: status
    tuple val(meta), path("streme.log"), emit: logs
    tuple val(meta), path("streme_versions.yml"), emit: versions

    script:
    stremeSampleId(meta)

    """
    set -euo pipefail
    if [[ ! -s "foreground.fa" ]]; then
        mkdir -p "streme"
        : > "streme/streme.txt"
        printf 'status\\tno_peaks\\n' > "streme_status.tsv"
        printf '%s\\n' 'no peaks; STREME skipped' > "streme.log"
    else
        if [[ ! -s "background.fa" ]]; then
            printf '%s\\n' \
                'foreground peaks exist but shuffled background FASTA is empty' >&2
            exit 2
        fi
        streme --oc "streme" \
            --dna \
            --p "foreground.fa" \
            --n "background.fa" \
            > "streme.log" 2>&1
        printf 'status\\tcomputed\\n' > "streme_status.tsv"
    fi

    printf 'STREME:\\n  streme: ' > "streme_versions.yml"
    streme --version 2>&1 | sed -n '1p' >> "streme_versions.yml"
    """
}
