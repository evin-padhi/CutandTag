process COUNT_DIFFERENTIAL_FRAGMENTS {
    tag "${assay_target} fragment counts"
    label 'process_heavy'

    conda "${projectDir}/envs/subread.yml"
    container 'quay.io/biocontainers/subread:2.1.1--h577a1d6_0'

    publishDir path: { "${params.outdir}/differential_binding/${assay_target}" },
        mode: 'copy',
        overwrite: true

    input:
    tuple val(assay_target), val(sample_ids),
        path(bam_files, stageAs: 'bams??/*', arity: '1..*'),
        path(consensus_bed, stageAs: 'consensus_peaks.bed'),
        path(consensus_saf, stageAs: 'consensus_peaks.saf')

    output:
    tuple val(assay_target), path("fragment_counts.tsv"), emit: counts
    tuple val(assay_target), path("count_differential_fragments.log"), emit: logs
    tuple val(assay_target), path("count_differential_fragments_versions.yml"),
        emit: versions

    script:
    def sampleIds = sample_ids.join(' ')
    def quoteShell = { value ->
        "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }
    def bamSourceArgs = bam_files.collect(quoteShell).join(' ')

    """
    set -euo pipefail
    export LC_ALL=C
    sample_ids=( ${sampleIds} )
    bam_sources=( ${bamSourceArgs} )
    if [[ \${#bam_sources[@]} -ne \${#sample_ids[@]} ]]; then
        printf 'Expected %s target BAMs; found %s\\n' \\
            "\${#sample_ids[@]}" "\${#bam_sources[@]}" \\
            | tee "count_differential_fragments.log" >&2
        exit 1
    fi

    mkdir -p count_bams
    for index in "\${!sample_ids[@]}"; do
        ln -s "\$PWD/\${bam_sources[index]}" \\
            "count_bams/\${sample_ids[index]}.bam"
    done
    bam_paths=(count_bams/*.bam)
    if [[ \${#bam_paths[@]} -ne \${#sample_ids[@]} ]]; then
        printf 'Could not stage one named BAM per target sample\\n' \\
            | tee "count_differential_fragments.log" >&2
        exit 1
    fi

    featureCounts \\
        -T "${task.cpus}" \\
        -p \\
        --countReadPairs \\
        -B \\
        -C \\
        -s 0 \\
        -F SAF \\
        -a "consensus_peaks.saf" \\
        -o "featurecounts_raw.tsv" \\
        "\${bam_paths[@]}" \\
        > "count_differential_fragments.log" 2>&1

    awk -v sample_ids="${sampleIds}" '
        BEGIN {
            FS = OFS = "\\t"
            sample_count = split(sample_ids, ids, " ")
            printf "peak_id\\tchrom\\tstart\\tend"
            for (i = 1; i <= sample_count; i++) printf "\\t%s", ids[i]
            printf "\\n"
        }
        /^#/ || \$1 == "Geneid" { next }
        NF >= 6 + sample_count {
            printf "%s\\t%s\\t%d\\t%d", \$1, \$2, \$3 - 1, \$4
            for (i = 1; i <= sample_count; i++) printf "\\t%s", \$(6 + i)
            printf "\\n"
        }
    ' "featurecounts_raw.tsv" > "fragment_counts.tsv"

    printf 'COUNT_DIFFERENTIAL_FRAGMENTS:\\n  subread: ' \\
        > "count_differential_fragments_versions.yml"
    featureCounts -v 2>&1 | sed -n '1p' \\
        >> "count_differential_fragments_versions.yml"
    """
}
