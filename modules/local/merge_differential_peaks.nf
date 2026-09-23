process MERGE_DIFFERENTIAL_PEAKS {
    tag "${assay_target} consensus regions"
    label 'process_standard'

    conda "${projectDir}/envs/bedtools.yml"
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'

    publishDir "${params.outdir}/differential_binding/${assay_target}",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename == 'consensus_peaks.saf' ? null : filename
        }

    input:
    tuple val(assay_target), val(sample_ids),
        path(peak_files, stageAs: 'peaks??/*', arity: '1..*')

    output:
    tuple val(assay_target), val(sample_ids),
        path("consensus_peaks.bed"), path("consensus_peaks.saf"),
        emit: consensus
    tuple val(assay_target), path("merge_differential_peaks.log"), emit: logs
    tuple val(assay_target), path("merge_differential_peaks_versions.yml"),
        emit: versions

    script:
    def sampleIds = sample_ids.join(' ')

    """
    set -euo pipefail
    export LC_ALL=C
    sample_ids=( ${sampleIds} )
    mapfile -t peak_sources < <(
        find peaks?? -maxdepth 2 -type f -name 'final.broadPeak' -print | sort
    )
    if [[ \${#peak_sources[@]} -ne \${#sample_ids[@]} ]]; then
        printf 'Expected %s target peak files; found %s\\n' \\
            "\${#sample_ids[@]}" "\${#peak_sources[@]}" \\
            | tee "merge_differential_peaks.log" >&2
        exit 1
    fi

    : > all_target_peaks.bed
    for index in "\${!sample_ids[@]}"; do
        awk -F '\\t' 'BEGIN { OFS = "\\t" }
            NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 < $3 {
                print $1, $2, $3
            }' "\${peak_sources[$index]}" >> all_target_peaks.bed
    done
    if [[ ! -s all_target_peaks.bed ]]; then
        printf 'No target peak intervals were called across assay %s\\n' \\
            "${assay_target}" | tee "merge_differential_peaks.log" >&2
        exit 1
    fi

    bedtools sort -i all_target_peaks.bed \\
        | bedtools merge -i - > "consensus_peaks.bed" \\
        2> "merge_differential_peaks.log"
    if [[ ! -s consensus_peaks.bed ]]; then
        printf 'Consensus peak set is empty for assay %s\\n' \\
            "${assay_target}" | tee -a "merge_differential_peaks.log" >&2
        exit 1
    fi

    awk 'BEGIN { OFS = "\\t"; print "GeneID", "Chr", "Start", "End", "Strand" }
        { printf "peak_%07d\\t%s\\t%d\\t%d\\t+\\n", NR, $1, $2 + 1, $3 }' \\
        consensus_peaks.bed > "consensus_peaks.saf"

    printf 'MERGE_DIFFERENTIAL_PEAKS:\\n  bedtools: ' \\
        > "merge_differential_peaks_versions.yml"
    bedtools --version 2>&1 | sed -n '1p' \\
        >> "merge_differential_peaks_versions.yml"
    """
}
