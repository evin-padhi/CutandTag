process SAMTOOLS_METRICS {
    tag "${meta.sample_id}"
    label 'process_standard'

    conda "${projectDir}/envs/samtools.yml"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    publishDir "${params.outdir}/qc/library/${meta.sample_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(analysis_bam, stageAs: 'metrics.analysis.bam'),
        path(analysis_bai, stageAs: 'metrics.analysis.bam.bai')

    output:
    tuple val(meta),
        path("flagstat.txt"),
        path("alignment.stats.txt"),
        path("idxstats.tsv"),
        path("insert_size.tsv"),
        path("duplicate_metrics.json"),
        emit: metrics
    tuple val(meta), path("samtools_metrics_versions.yml"), emit: versions

    script:
    """
    set -euo pipefail
    samtools flagstat \
        -@ "${task.cpus}" \
        "metrics.analysis.bam" \
        > "flagstat.txt"
    samtools stats \
        "metrics.analysis.bam" \
        > "alignment.stats.txt"
    samtools idxstats \
        "metrics.analysis.bam" \
        > "idxstats.tsv"

    printf 'record_type\\tinsert_size\\ttotal_pairs\\tinward_pairs\\toutward_pairs\\tother_pairs\\n' \
        > "insert_size.tsv"
    awk -F '\\t' '\$1 == "IS" { print \$0 }' \
        "alignment.stats.txt" \
        >> "insert_size.tsv"

    samtools collate \
        -@ "${task.cpus}" \
        -u \
        -O \
        "metrics.analysis.bam" \
        | samtools fixmate -m -@ "${task.cpus}" -u - - \
        | samtools sort -@ "${task.cpus}" -u -o "markdup.input.bam" -
    samtools markdup \
        --json \
        -s \
        -f "duplicate_metrics.json" \
        -@ "${task.cpus}" \
        -O bam \
        "markdup.input.bam" \
        "/dev/null"

    printf 'SAMTOOLS_METRICS:\\n  samtools: ' > "samtools_metrics_versions.yml"
    samtools --version 2>&1 | sed -n '1p' >> "samtools_metrics_versions.yml"
    """
}
