process FILTER_BLACKLIST {
    tag "${meta.sample_id}"
    label 'process_standard'

    conda "${projectDir}/envs/bedtools.yml"
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'

    publishDir "${params.outdir}/peaks/${meta.sample_id}/broad/final",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta), path(raw_broad_peak, stageAs: 'raw.broadPeak')
    path(blacklist_files, stageAs: 'blacklist/regions*.bed', arity: '0..1')

    output:
    tuple val(meta), path("final.broadPeak"), emit: peaks
    tuple val(meta), path("filter_blacklist.log"), emit: logs
    tuple val(meta), path("filter_blacklist_versions.yml"), emit: versions

    script:
    def suppliedBlacklists
    if (blacklist_files == null) {
        suppliedBlacklists = []
    } else if (blacklist_files instanceof java.util.Collection) {
        suppliedBlacklists = blacklist_files as List
    } else {
        suppliedBlacklists = [blacklist_files]
    }
    if (suppliedBlacklists.size() > 1) {
        throw new IllegalArgumentException(
            "FILTER_BLACKLIST accepts at most one blacklist, got " +
            suppliedBlacklists.size()
        )
    }
    def hasBlacklist = suppliedBlacklists.size() == 1
    def finalizeCommand = hasBlacklist
        ? '''
    bedtools intersect \
        -v \
        -a "raw.broadPeak" \
        -b "blacklist/regions.bed" \
        > "final.broadPeak" \
        2> "filter_blacklist.log"
    '''
        : '''
    cp "raw.broadPeak" "final.broadPeak"
    printf 'blacklist absent; raw broadPeak retained unchanged\\n' \
        > "filter_blacklist.log"
    '''
    def versionCommand = hasBlacklist
        ? '''
    printf 'FILTER_BLACKLIST:\\n  bedtools: ' \
        > "filter_blacklist_versions.yml"
    bedtools --version 2>&1 | sed -n '1p' \
        >> "filter_blacklist_versions.yml"
    '''
        : '''
    printf 'FILTER_BLACKLIST:\\n  bedtools: not run (blacklist absent)\\n' \
        > "filter_blacklist_versions.yml"
    '''

    """
    set -euo pipefail
    ${finalizeCommand}
    ${versionCommand}
    """
}
