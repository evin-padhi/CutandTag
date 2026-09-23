include { MERGE_DIFFERENTIAL_PEAKS } from '../../modules/local/merge_differential_peaks'
include { COUNT_DIFFERENTIAL_FRAGMENTS } from '../../modules/local/count_differential_fragments'

workflow DIFFERENTIAL_BINDING {
    take:
    filtered_bams
    narrow_peaks
    blacklist

    main:
    target_bams = filtered_bams
        .filter { meta, bam, bai -> meta.is_control == false }
        .map { meta, bam, bai ->
            tuple(meta.sample_id.toString(), meta.assay_target.toString(), bam)
        }

    target_peaks = narrow_peaks.map { meta, peaks ->
        tuple(meta.sample_id.toString(), meta.assay_target.toString(),
            meta.condition.toString(), peaks)
    }

    paired_targets = target_bams.join(target_peaks).map {
            sample_id, bam_assay, bam, peak_assay, condition, peaks ->
            if (bam_assay != peak_assay) {
                throw new IllegalStateException(
                    "assay target mismatch for ${sample_id}: " +
                    "BAM is ${bam_assay}; peak file is ${peak_assay}"
                )
            }
            tuple(bam_assay, sample_id, condition, bam, peaks)
    }

    assay_groups = paired_targets
        .groupTuple(by: 0)
        .map { assay_target, raw_sample_ids, raw_conditions, raw_bams, raw_peaks ->
            def rows = (0..<raw_sample_ids.size()).collect { index ->
                [
                    sample_id: raw_sample_ids[index].toString(),
                    condition: raw_conditions[index].toString(),
                    bam: raw_bams[index],
                    peaks: raw_peaks[index],
                ]
            }.sort { left, right -> left.sample_id <=> right.sample_id }
            tuple(
                assay_target.toString(),
                rows.collect { row -> row.sample_id },
                rows.collect { row -> row.condition },
                rows.collect { row -> row.bam },
                rows.collect { row -> row.peaks }
            )
        }

    reusable_blacklist = blacklist
        .reduce([files: []]) { holder, path -> [files: holder.files + [path]] }
        .map { holder -> holder.files }
    peak_inputs = assay_groups
        .combine(reusable_blacklist)
        .map { assay_target, sample_ids, conditions, bams, peaks, blacklists ->
            tuple(assay_target, sample_ids, conditions, peaks, blacklists)
        }
    bam_inputs = assay_groups.map { assay_target, sample_ids, conditions, bams, peaks ->
        tuple(assay_target, sample_ids, bams)
    }

    MERGE_DIFFERENTIAL_PEAKS(peak_inputs)

    count_inputs = bam_inputs
        .join(MERGE_DIFFERENTIAL_PEAKS.out.consensus, by: 0)
        .map { assay_target, sample_ids, bams, peak_sample_ids, bed, saf ->
            if (sample_ids != peak_sample_ids) {
                throw new IllegalStateException(
                    "sample ordering changed while preparing ${assay_target} " +
                    "differential-binding inputs"
                )
            }
            tuple(assay_target, sample_ids, bams, bed, saf)
        }
    COUNT_DIFFERENTIAL_FRAGMENTS(count_inputs)

    emit:
    consensus_peaks = MERGE_DIFFERENTIAL_PEAKS.out.consensus.map {
        assay_target, sample_ids, bed, saf -> tuple(assay_target, bed)
    }
    fragment_counts = COUNT_DIFFERENTIAL_FRAGMENTS.out.counts
    logs = MERGE_DIFFERENTIAL_PEAKS.out.logs.mix(
        COUNT_DIFFERENTIAL_FRAGMENTS.out.logs
    )
    versions = MERGE_DIFFERENTIAL_PEAKS.out.versions.mix(
        COUNT_DIFFERENTIAL_FRAGMENTS.out.versions
    )
}
