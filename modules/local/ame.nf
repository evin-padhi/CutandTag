def ameSampleId(meta) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    def sampleId = meta instanceof Map ? meta.sample_id?.toString() : null
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "AME sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    sampleId
}

process AME {
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
    tuple val(meta), path("ame"), emit: results
    tuple val(meta), path("ame_status.tsv"), emit: status
    tuple val(meta), path("ame.log"), emit: logs
    tuple val(meta), path("ame_versions.yml"), emit: versions

    script:
    ameSampleId(meta)

    """
    set -euo pipefail
    if [[ ! -s "foreground.fa" ]]; then
        mkdir -p "ame"
        printf 'rank\\tmotif_ID\\tmotif_Alt_ID\\tp-value\\tE-value\\tpos\\tneg\\n' \
            > "ame/ame.tsv"
        printf '1\\t__NO_PEAKS__\\t__NO_PEAKS__\\t1\\t1\\t0\\t0\\n' \
            >> "ame/ame.tsv"
        printf 'status\\tno_peaks\\n' > "ame_status.tsv"
        printf '%s\\n' 'no peaks; AME skipped' > "ame.log"
    else
        if [[ ! -s "background.fa" ]]; then
            printf '%s\\n' \
                'foreground peaks exist but shuffled background FASTA is empty' >&2
            exit 2
        fi
        mkdir -p "ame"
        ame --text \
            --control "background.fa" \
            --method fisher \
            --scoring avg \
            --evalue-report-threshold 1e300 \
            "foreground.fa" \
            "motifs.meme" \
            > "ame/ame.tsv" \
            2> "ame.log"
        if [[ ! -f "ame/ame.tsv" ]]; then
            printf '%s\\n' 'AME completed without ame/ame.tsv' >&2
            exit 2
        fi
        {
            printf '%s\\n' '# motif_qc_complete_database=true'
            cat "ame/ame.tsv"
        } > "ame/ame.complete.tsv"
        mv "ame/ame.complete.tsv" "ame/ame.tsv"
        printf 'status\\tcomputed\\n' > "ame_status.tsv"
    fi

    printf 'AME:\\n  ame: ' > "ame_versions.yml"
    ame --version 2>&1 | sed -n '1p' >> "ame_versions.yml"
    """
}
