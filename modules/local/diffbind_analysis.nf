process DIFFBIND_ANALYSIS {
    tag "${assay_target} differential binding"
    label 'process_heavy'
    errorStrategy 'ignore'

    conda "${projectDir}/envs/diffbind.yml"
    container 'ghcr.io/evin-padhi/cutandtag-diffbind:diffbind-0.1.0'

    publishDir path: { "${params.outdir}/differential_binding/${assay_target}" },
        mode: 'copy',
        overwrite: true

    input:
    tuple val(assay_target), val(sample_rows),
        path(bam_files, stageAs: 'bams??/filtered.bam', arity: '1..*'),
        path(bai_files, stageAs: 'bams??/filtered.bam.bai', arity: '1..*'),
        path(consensus_bed, stageAs: 'consensus_peaks.bed'),
        path(analysis_script, stageAs: 'diffbind_analysis.R')

    output:
    tuple val(assay_target), path('diffbind_results.tsv'), emit: results
    tuple val(assay_target), path('diffbind_comparison_summary.tsv'), emit: summary
    tuple val(assay_target), path('diffbind_mode_correlations.tsv'), emit: mode_correlations
    tuple val(assay_target), path('replicate_correlation.tsv'), emit: correlations
    tuple val(assay_target), path('consensus_peak_ids.tsv'), emit: peak_ids
    tuple val(assay_target), path('*.pdf', arity: '1..*'), emit: plots
    tuple val(assay_target), path('diffbind_analysis.log'), emit: logs
    tuple val(assay_target), path('diffbind_analysis_versions.yml'), emit: versions

    script:
    def targetRows = sample_rows.findAll { row -> !row.is_control }
    if (targetRows.isEmpty()) {
        throw new IllegalArgumentException(
            "DIFFBIND_ANALYSIS received no target samples for ${assay_target}"
        )
    }
    def quoteShell = { value ->
        "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }
    def sheetCommands = targetRows.collect { row ->
        def control = sample_rows.find { candidate ->
            candidate.sample_id == row.control_id
        }
        if (control == null || !control.is_control) {
            throw new IllegalArgumentException(
                "No matched IgG row for target ${row.sample_id} (${row.control_id})"
            )
        }
        def targetBam = bam_files[row.bam_index].toString()
        def controlBam = bam_files[control.bam_index].toString()
        def fields = [row.sample_id, row.condition, targetBam, row.control_id, controlBam]
        "printf '%s\\t%s\\t%s\\t%s\\t%s\\n' ${fields.collect(quoteShell).join(' ')} >> sample_sheet.tsv"
    }.join('\n')

    """
    set -euo pipefail
    export LC_ALL=C
    printf '%s\\n' 'sample_id	condition	bam_reads	control_id	bam_control' \
        > sample_sheet.tsv
    ${sheetCommands}

    Rscript "diffbind_analysis.R" \
        --sample-sheet sample_sheet.tsv \
        --consensus-bed "consensus_peaks.bed" \
        --assay "${assay_target}" \
        --fdr 0.05 \
        --outdir results

    cp results/diffbind_results.tsv .
    cp results/diffbind_comparison_summary.tsv .
    cp results/diffbind_mode_correlations.tsv .
    cp results/replicate_correlation.tsv .
    cp results/consensus_peak_ids.tsv .
    cp results/diffbind_analysis.log .
    cp results/diffbind_analysis_versions.yml .
    cp results/*.pdf .
    """
}
