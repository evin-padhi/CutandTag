include { BAM_TO_FRAGMENTS } from '../../modules/local/bam_to_fragments'
include { FILTERED_BAM_QC } from '../../modules/local/filtered_bam_qc'
include { PEAK_QC } from '../../modules/local/peak_qc'
include { TSS_ENRICHMENT } from '../../modules/local/tss_enrichment'
include { DEMUX_QC_CUSTOM; LIBRARY_QC_CUSTOM; MOTIF_QC_CUSTOM; MULTIQC } from '../../modules/local/multiqc'

def validateQcMeta(rawMeta, context) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    if (!(rawMeta instanceof Map)) {
        throw new IllegalArgumentException(
            "${context} metadata must be a map"
        )
    }
    def meta = new LinkedHashMap(rawMeta)
    def sampleId = meta.sample_id?.toString()
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "${context} sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    if (!(meta.is_control instanceof Boolean)) {
        throw new IllegalArgumentException(
            "${context} is_control must be a boolean for sample ${sampleId}"
        )
    }
    meta.sample_id = sampleId
    meta
}

def singletonOptional(rows, label) {
    if (rows.size() > 1) {
        throw new IllegalArgumentException(
            "QC accepts at most one ${label}, got ${rows.size()}"
        )
    }
    rows
}

workflow QC {
    take:
    filtered_bams
    final_broad_peaks
    coverage
    library_metrics
    demultiplex_metrics
    fastqc_reports
    motif_metrics
    gtf
    tss_bed

    main:
    safe_filtered_bams = filtered_bams.map { meta, bam, bai ->
        tuple(validateQcMeta(meta, 'filtered BAM'), bam, bai)
    }
    safe_final_broad_peaks = final_broad_peaks.map { meta, peaks ->
        def safeMeta = validateQcMeta(meta, 'final broad peak')
        if (safeMeta.is_control) {
            throw new IllegalArgumentException(
                "IgG control ${safeMeta.sample_id} cannot enter default FRiP"
            )
        }
        tuple(safeMeta, peaks)
    }
    safe_coverage = coverage.map { meta, bigwig ->
        tuple(validateQcMeta(meta, 'coverage'), bigwig)
    }
    safe_library_metrics = library_metrics.map {
        meta, flagstat, stats, idxstats, insertSize, duplicateMetrics ->
        tuple(
            validateQcMeta(meta, 'library QC'),
            flagstat,
            stats,
            idxstats,
            insertSize,
            duplicateMetrics
        )
    }
    /*
     * Task 10 contract: tuple(meta, expected_motif_qc_tsv).
     * sample_id comes only from metadata; motif_qc.py's metric/value TSV does
     * not duplicate pipeline identity inside its content.
     */
    safe_motif_metrics = motif_metrics.map { meta, motifTsv ->
        def safeMeta = validateQcMeta(meta, 'motif QC')
        if (safeMeta.is_control) {
            throw new IllegalArgumentException(
                "IgG control ${safeMeta.sample_id} cannot have expected-motif QC"
            )
        }
        tuple(safeMeta, motifTsv)
    }

    /*
     * Only non-control filtered BAMs can enter fragment conversion and FRiP.
     * Keyed fail-on-mismatch joining makes a missing or extraneous final
     * broadPeak an error instead of silently dropping a target.
     */
    target_filtered_bams = safe_filtered_bams.filter {
        meta, bam, bai -> !meta.is_control
    }
    keyed_target_bams = target_filtered_bams.map { meta, bam, bai ->
        tuple(meta.sample_id, meta, bam, bai)
    }
    keyed_final_broad_peaks = safe_final_broad_peaks.map { meta, peaks ->
        tuple(meta.sample_id, meta, peaks)
    }
    matched_targets = keyed_target_bams.join(
        keyed_final_broad_peaks,
        by: 0,
        failOnDuplicate: true,
        failOnMismatch: true
    ).map {
        sampleId, bamMeta, bam, bai, peakMeta, peaks ->
        if (bamMeta.sample_id != peakMeta.sample_id) {
            throw new IllegalStateException(
                "QC target metadata mismatch for ${sampleId}"
            )
        }
        tuple(bamMeta, bam, bai, peaks)
    }

    fragment_inputs = matched_targets.map { meta, bam, bai, peaks ->
        tuple(meta, bam, bai)
    }
    matched_peak_inputs = matched_targets.map { meta, bam, bai, peaks ->
        tuple(meta.sample_id, meta, peaks)
    }

    BAM_TO_FRAGMENTS(fragment_inputs)

    keyed_fragments = BAM_TO_FRAGMENTS.out.fragments.map {
        meta, fragments -> tuple(meta.sample_id, meta, fragments)
    }
    peak_qc_inputs = keyed_fragments.join(
        matched_peak_inputs,
        by: 0,
        failOnDuplicate: true,
        failOnMismatch: true
    ).map {
        sampleId, fragmentMeta, fragments, peakMeta, peaks ->
        tuple(fragmentMeta, fragments, peaks)
    }
    PEAK_QC(peak_qc_inputs)

    /*
     * Count reads/fragments in each filtered BAM with SAMtools, then
     * fail-on-mismatch join those counts to Task 7 library metrics. This is
     * the numerator for MAPQ-filtered fraction and cannot be inferred from
     * raw flagstat/stats files.
     */
    FILTERED_BAM_QC(safe_filtered_bams)
    keyed_library_metrics = safe_library_metrics.map {
        meta, flagstat, stats, idxstats, insertSize, duplicateMetrics ->
        tuple(
            meta.sample_id,
            meta,
            flagstat,
            stats,
            idxstats,
            insertSize,
            duplicateMetrics
        )
    }
    keyed_filtered_bam_qc = FILTERED_BAM_QC.out.metrics.map {
        meta, filteredMetrics ->
        tuple(meta.sample_id, meta, filteredMetrics)
    }
    library_qc_inputs = keyed_library_metrics.join(
        keyed_filtered_bam_qc,
        by: 0,
        failOnDuplicate: true,
        failOnMismatch: true
    ).map {
        sampleId,
        metricMeta,
        flagstat,
        stats,
        idxstats,
        insertSize,
        duplicateMetrics,
        filteredMeta,
        filteredMetrics ->
        if (
            metricMeta.sample_id != filteredMeta.sample_id ||
            metricMeta.is_control != filteredMeta.is_control
        ) {
            throw new IllegalStateException(
                "library QC metadata mismatch for ${sampleId}"
            )
        }
        tuple(
            metricMeta,
            flagstat,
            stats,
            idxstats,
            insertSize,
            duplicateMetrics,
            filteredMetrics
        )
    }

    /*
     * Each optional path queue is collected into a value, including when it
     * is empty. A supplied TSS BED takes precedence over GTF as documented.
     * Therefore absent annotations close cleanly and never deadlock coverage.
     */
    reusable_gtf = gtf.collect(flat: false)
        .ifEmpty { ignored -> [] }
        .map { rows -> [files: singletonOptional(rows, 'GTF')] }
    reusable_tss_bed = tss_bed.collect(flat: false)
        .ifEmpty { ignored -> [] }
        .map { rows -> [files: singletonOptional(rows, 'TSS BED')] }
    annotation_choice = reusable_gtf.combine(reusable_tss_bed).map {
        gtf_holder, tss_holder ->
        def gtf_rows = gtf_holder.files
        def tss_rows = tss_holder.files
        def annotation_mode = tss_rows ? 'bed' : (gtf_rows ? 'gtf' : 'none')
        def annotation_files = tss_rows ?: gtf_rows
        tuple(annotation_mode, annotation_files)
    }
    annotation_status = annotation_choice.map { annotation_mode, files ->
        annotation_mode == 'none'
            ? 'skipped_no_annotation'
            : "computed_${annotation_mode}"
    }

    tss_inputs = safe_coverage.combine(annotation_choice)
        .filter {
            meta, bigwig, annotation_mode, annotation_files ->
            annotation_mode != 'none'
        }
        .map {
            meta, bigwig, annotation_mode, annotation_files ->
            tuple(meta, bigwig, annotation_mode, annotation_files)
        }
    TSS_ENRICHMENT(tss_inputs)

    /*
     * IgG and target libraries both receive these custom library tables.
     * Demultiplex tables are physical-library keyed. Peak tables remain
     * target-only because PEAK_QC is fed by target_filtered_bams above.
     */
    DEMUX_QC_CUSTOM(demultiplex_metrics)
    LIBRARY_QC_CUSTOM(library_qc_inputs)
    MOTIF_QC_CUSTOM(safe_motif_metrics)

    fastqc_files = fastqc_reports.flatMap { meta, html, zip ->
        def files = []
        [html, zip].each { value ->
            if (value instanceof java.util.Collection) {
                files.addAll(value)
            } else if (value != null) {
                files.add(value)
            }
        }
        files
    }.reduce([files: []]) { holder, path ->
        [files: holder.files + [path]]
    }.map { holder -> holder.files }
    demux_custom_files = DEMUX_QC_CUSTOM.out.custom
        .map { meta, custom -> custom }
        .reduce([files: []]) { holder, path ->
            [files: holder.files + [path]]
        }
        .map { holder -> holder.files }
    library_custom_files = LIBRARY_QC_CUSTOM.out.custom
        .map { meta, custom, insertDistribution -> custom }
        .reduce([files: []]) { holder, path ->
            [files: holder.files + [path]]
        }
        .map { holder -> holder.files }
    insert_size_files = LIBRARY_QC_CUSTOM.out.custom
        .map { meta, custom, insertDistribution -> insertDistribution }
        .reduce([files: []]) { holder, path ->
            [files: holder.files + [path]]
        }
        .map { holder -> holder.files }
    peak_qc_files = PEAK_QC.out.qc
        .map { meta, json, tsv, histogram, perPeak -> tsv }
        .reduce([files: []]) { holder, path ->
            [files: holder.files + [path]]
        }
        .map { holder -> holder.files }
    motif_metric_files = MOTIF_QC_CUSTOM.out.custom
        .map { meta, custom -> custom }
        .reduce([files: []]) { holder, path ->
            [files: holder.files + [path]]
        }
        .map { holder -> holder.files }
    tss_status_files = TSS_ENRICHMENT.out.profiles
        .map { meta, bed, matrix, matrixTable, profile, status -> status }
        .reduce([files: []]) { holder, path ->
            [files: holder.files + [path]]
        }
        .map { holder -> holder.files }

    MULTIQC(
        fastqc_files,
        demux_custom_files,
        library_custom_files,
        insert_size_files,
        peak_qc_files,
        motif_metric_files,
        tss_status_files,
        annotation_status
    )

    versions_ch = BAM_TO_FRAGMENTS.out.versions.mix(
        FILTERED_BAM_QC.out.versions,
        PEAK_QC.out.versions,
        TSS_ENRICHMENT.out.versions,
        MULTIQC.out.versions
    )

    emit:
    target_fragments = BAM_TO_FRAGMENTS.out.fragments
    target_qc = PEAK_QC.out.qc
    combined_summary = MULTIQC.out.combined_summary
    tss_profiles = TSS_ENRICHMENT.out.profiles
    multiqc_report = MULTIQC.out.report
    multiqc_data = MULTIQC.out.data
    multiqc_custom_content = MULTIQC.out.custom_content
    versions = versions_ch
}
