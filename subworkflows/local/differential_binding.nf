include { MERGE_DIFFERENTIAL_PEAKS } from '../../modules/local/merge_differential_peaks'
include { COUNT_DIFFERENTIAL_FRAGMENTS } from '../../modules/local/count_differential_fragments'
include { DIFFBIND_ANALYSIS } from '../../modules/local/diffbind_analysis'

workflow DIFFERENTIAL_BINDING {
    take:
    filtered_bams
    narrow_peaks
    blacklist

    main:
    all_bams = filtered_bams
        .map { meta, bam, bai ->
            tuple(new LinkedHashMap(meta), bam, bai)
        }
        .collect()

    target_bams = filtered_bams
        .filter { meta, bam, bai -> !meta.is_control }
        .map { meta, bam, bai ->
            tuple(meta.sample_id.toString(), new LinkedHashMap(meta), bam, bai)
        }

    target_peaks = narrow_peaks.map { meta, peaks ->
        tuple(meta.sample_id.toString(), meta.assay_target.toString(),
            meta.condition.toString(), peaks)
    }

    paired_targets = target_bams.join(target_peaks).map {
            sample_id, meta, bam, bai, peak_assay, condition, peaks ->
            if (meta.assay_target.toString() != peak_assay) {
                throw new IllegalStateException(
                    "assay target mismatch for ${sample_id}: " +
                    "BAM is ${meta.assay_target}; peak file is ${peak_assay}"
                )
            }
            tuple(meta.assay_target.toString(), sample_id, condition, meta, bam, bai, peaks)
    }

    assay_groups = paired_targets
        .groupTuple(by: 0)
        .map { assay_target, raw_sample_ids, raw_conditions, raw_metas,
               raw_bams, raw_bais, raw_peaks ->
            def rows = (0..<raw_sample_ids.size()).collect { index ->
                [
                    sample_id: raw_sample_ids[index].toString(),
                    condition: raw_conditions[index].toString(),
                    control_id: raw_metas[index].control_id.toString(),
                    meta: raw_metas[index],
                    bam: raw_bams[index],
                    bai: raw_bais[index],
                    peaks: raw_peaks[index],
                ]
            }.sort { left, right -> left.sample_id <=> right.sample_id }
            tuple(assay_target.toString(), rows)
        }

    reusable_blacklist = blacklist
        .reduce([files: []]) { holder, path -> [files: holder.files + [path]] }
        .map { holder -> holder.files }
    peak_inputs = assay_groups
        .combine(reusable_blacklist)
        .map { assay_target, rows, blacklists ->
            tuple(
                assay_target,
                rows.collect { row -> row.sample_id },
                rows.collect { row -> row.condition },
                rows.collect { row -> row.peaks },
                blacklists
            )
        }

    MERGE_DIFFERENTIAL_PEAKS(peak_inputs)

    count_inputs = assay_groups
        .join(MERGE_DIFFERENTIAL_PEAKS.out.consensus, by: 0)
        .map { assay_target, rows, peak_sample_ids, bed, saf ->
            def sample_ids = rows.collect { row -> row.sample_id }
            if (sample_ids != peak_sample_ids) {
                throw new IllegalStateException(
                    "sample ordering changed while preparing ${assay_target} " +
                    "differential-binding inputs"
                )
            }
            tuple(
                assay_target,
                sample_ids,
                rows.collect { row -> row.bam },
                bed,
                saf
            )
        }
    COUNT_DIFFERENTIAL_FRAGMENTS(count_inputs)

    all_bam_inventory = all_bams.map { records ->
        records.collectEntries { meta, bam, bai ->
            [(meta.sample_id.toString()): [meta: meta, bam: bam, bai: bai]]
        }
    }
    diffbind_inputs = assay_groups
        .combine(all_bam_inventory)
        .join(MERGE_DIFFERENTIAL_PEAKS.out.consensus, by: 0)
        .map { assay_target, targets, inventory,
               peak_sample_ids, consensus_bed, consensus_saf ->
            def target_ids = targets.collect { row -> row.sample_id }
            if (target_ids != peak_sample_ids) {
                throw new IllegalStateException(
                    "target sample ordering changed while preparing ${assay_target} " +
                    "DiffBind inputs"
                )
            }
            def target_controls = targets.collect { row -> row.control_id }.unique()
            def control_rows = target_controls.collect { control_id ->
                def control = inventory[control_id]
                if (control == null || !control.meta.is_control ||
                    control.meta.assay_target.toString() != 'IgG') {
                    throw new IllegalStateException(
                        "target control ${control_id} is missing or is not an IgG sample " +
                        "for assay ${assay_target}"
                    )
                }
                [
                    sample_id: control_id,
                    condition: '',
                    control_id: '',
                    is_control: true,
                    bam: control.bam,
                    bai: control.bai
                ]
            }
            def all_rows = (targets.collect { row ->
                [
                    sample_id: row.sample_id,
                    condition: row.condition,
                    control_id: row.control_id,
                    is_control: false,
                    bam: row.bam,
                    bai: row.bai
                ]
            } + control_rows).unique { row -> row.sample_id }
                .sort { left, right -> left.sample_id <=> right.sample_id }
            def sample_rows = all_rows.withIndex().collect { row, index ->
                new LinkedHashMap(row) + [bam_index: index]
            }
            tuple(
                assay_target,
                sample_rows,
                all_rows.collect { row -> row.bam },
                all_rows.collect { row -> row.bai },
                consensus_bed
            )
        }
    analysis_script = Channel.fromPath(
        "${projectDir}/bin/diffbind_analysis.R",
        checkIfExists: true
    )
    DIFFBIND_ANALYSIS(diffbind_inputs.combine(analysis_script))

    emit:
    consensus_peaks = MERGE_DIFFERENTIAL_PEAKS.out.consensus.map {
        assay_target, sample_ids, bed, saf -> tuple(assay_target, bed)
    }
    fragment_counts = COUNT_DIFFERENTIAL_FRAGMENTS.out.counts
    diffbind_results = DIFFBIND_ANALYSIS.out.results
    diffbind_summary = DIFFBIND_ANALYSIS.out.summary
    replicate_correlations = DIFFBIND_ANALYSIS.out.correlations
    diffbind_peak_ids = DIFFBIND_ANALYSIS.out.peak_ids
    diffbind_plots = DIFFBIND_ANALYSIS.out.plots
    logs = MERGE_DIFFERENTIAL_PEAKS.out.logs.mix(
        COUNT_DIFFERENTIAL_FRAGMENTS.out.logs,
        DIFFBIND_ANALYSIS.out.logs
    )
    versions = MERGE_DIFFERENTIAL_PEAKS.out.versions.mix(
        COUNT_DIFFERENTIAL_FRAGMENTS.out.versions,
        DIFFBIND_ANALYSIS.out.versions
    )
}
