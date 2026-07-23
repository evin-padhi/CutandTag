def fimoSampleId(meta) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    def sampleId = meta instanceof Map ? meta.sample_id?.toString() : null
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "FIMO sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    sampleId
}

process FIMO {
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
    path motif_db, stageAs: 'motifs.meme'

    output:
    tuple val(meta), path("fimo"), emit: results
    tuple val(meta), path("fimo_status.tsv"), emit: status
    tuple val(meta), path("fimo.log"), emit: logs
    tuple val(meta), path("fimo_versions.yml"), emit: versions

    script:
    fimoSampleId(meta)

    """
    set -euo pipefail
    if [[ ! -s "foreground.fa" ]]; then
        mkdir -p "fimo"
        printf '%s\\n' \
            '# motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence' \
            > "fimo/fimo.tsv"
        printf 'status\\tno_peaks\\n' > "fimo_status.tsv"
        printf '%s\\n' 'no peaks; FIMO skipped' > "fimo.log"
    else
        fimo --oc "fimo" \
            "motifs.meme" "foreground.fa" \
            > "fimo.log" 2>&1
        if [[ ! -f "fimo/fimo.tsv" ]]; then
            printf '%s\\n' 'FIMO completed without fimo/fimo.tsv' >&2
            exit 2
        fi
        printf 'status\\tcomputed\\n' > "fimo_status.tsv"
    fi

    printf 'FIMO:\\n  fimo: ' > "fimo_versions.yml"
    fimo --version 2>&1 | sed -n '1p' >> "fimo_versions.yml"
    """
}
