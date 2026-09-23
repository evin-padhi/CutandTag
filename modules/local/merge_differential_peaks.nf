process MERGE_DIFFERENTIAL_PEAKS {
    tag "${assay_target} consensus regions"
    label 'process_standard'

    conda "${projectDir}/envs/bedtools.yml"
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'

    publishDir path: { "${params.outdir}/differential_binding/${assay_target}" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename == 'consensus_peaks.saf' ? null : filename
        }

    input:
    tuple val(assay_target), val(sample_ids), val(conditions),
        path(peak_files, stageAs: 'peaks??/*', arity: '1..*'),
        path(blacklist_files, stageAs: 'blacklist/regions*.bed')

    output:
    tuple val(assay_target), val(sample_ids),
        path("consensus_peaks.bed"), path("consensus_peaks.saf"),
        emit: consensus
    tuple val(assay_target), path("merge_differential_peaks.log"), emit: logs
    tuple val(assay_target), path("merge_differential_peaks_versions.yml"),
        emit: versions

    script:
    def sampleIds = sample_ids.join(' ')
    def conditionIds = conditions.join(' ')
    def quoteShell = { value ->
        "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }
    def peakSourceArgs = peak_files.collect(quoteShell).join(' ')
    def suppliedBlacklists = blacklist_files instanceof java.util.Collection
        ? blacklist_files as List
        : blacklist_files == null ? [] : [blacklist_files]
    if (suppliedBlacklists.size() > 1) {
        throw new IllegalArgumentException(
            "MERGE_DIFFERENTIAL_PEAKS accepts at most one blacklist, got " +
            suppliedBlacklists.size()
        )
    }
    def filterBlacklist = suppliedBlacklists.size() == 1
        ? '''
    bedtools intersect -v -a "merged_unfiltered.bed" \\
        -b "blacklist/regions.bed" > "consensus_peaks.bed" \\
        2>> "merge_differential_peaks.log"
    '''
        : '''
    cp "merged_unfiltered.bed" "consensus_peaks.bed"
    '''

    """
    set -euo pipefail
    export LC_ALL=C
    sample_ids=( ${sampleIds} )
    conditions=( ${conditionIds} )
    if [[ \${#conditions[@]} -ne \${#sample_ids[@]} ]]; then
        printf 'Sample and condition lists have different lengths\\n' \\
            | tee "merge_differential_peaks.log" >&2
        exit 1
    fi
    peak_sources=( ${peakSourceArgs} )
    if [[ \${#peak_sources[@]} -ne \${#sample_ids[@]} ]]; then
        printf 'Expected %s target peak files; found %s\\n' \\
            "\${#sample_ids[@]}" "\${#peak_sources[@]}" \\
            | tee "merge_differential_peaks.log" >&2
        exit 1
    fi

    condition_names=()
    for condition in "\${conditions[@]}"; do
        found=0
        for known_condition in "\${condition_names[@]-}"; do
            if [[ "\${condition}" == "\${known_condition}" ]]; then
                found=1
                break
            fi
        done
        if [[ \${found} -eq 0 ]]; then condition_names+=("\${condition}"); fi
    done
    : > reproducible_condition_peaks.bed
    condition_number=0
    for condition in "\${condition_names[@]}"; do
        condition_inputs=()
        for index in "\${!conditions[@]}"; do
            if [[ "\${conditions[index]}" == "\${condition}" ]]; then
                sorted_peak="condition_\${condition_number}_sample_\${index}.bed"
                bedtools sort -i "\${peak_sources[index]}" > "\${sorted_peak}"
                condition_inputs+=("\${sorted_peak}")
            fi
        done
        support=1
        if [[ \${#condition_inputs[@]} -gt 1 ]]; then support=2; fi
        bedtools multiinter -i "\${condition_inputs[@]}" \\
            | awk -v required="\${support}" 'BEGIN { OFS = "\\t" }
                \$4 >= required { print \$1, \$2, \$3 }' \\
            >> reproducible_condition_peaks.bed
        ((condition_number += 1))
    done 2> "merge_differential_peaks.log"
    if [[ ! -s reproducible_condition_peaks.bed ]]; then
        printf 'No reproducible narrow peaks were called for assay %s\\n' \\
            "${assay_target}" | tee -a "merge_differential_peaks.log" >&2
        exit 1
    fi

    bedtools sort -i reproducible_condition_peaks.bed \\
        | bedtools merge -i - > "merged_unfiltered.bed" \\
        2>> "merge_differential_peaks.log"
    ${filterBlacklist}
    if [[ ! -s consensus_peaks.bed ]]; then
        printf 'Consensus peak set is empty for assay %s\\n' \\
            "${assay_target}" | tee -a "merge_differential_peaks.log" >&2
        exit 1
    fi

    awk 'BEGIN { OFS = "\\t"; print "GeneID", "Chr", "Start", "End", "Strand" }
        { printf "peak_%07d\\t%s\\t%d\\t%d\\t+\\n", NR, \$1, \$2 + 1, \$3 }' \\
        consensus_peaks.bed > "consensus_peaks.saf"

    printf 'MERGE_DIFFERENTIAL_PEAKS:\\n  bedtools: ' \\
        > "merge_differential_peaks_versions.yml"
    bedtools --version 2>&1 | sed -n '1p' \\
        >> "merge_differential_peaks_versions.yml"
    """
}
