include { MACS2_BROAD } from '../../modules/local/macs2_broad'
include { MACS2_NARROW } from '../../modules/local/macs2_narrow'
include { FILTER_BLACKLIST } from '../../modules/local/filter_blacklist'

def isControlMeta(meta) {
    if (!(meta instanceof Map) || !(meta.is_control instanceof Boolean)) {
        throw new IllegalArgumentException(
            "peak metadata is_control must be a boolean for sample " +
            "${meta instanceof Map ? meta.sample_id : '<unknown>'}"
        )
    }
    meta.is_control
}

def validatePeakMeta(rawMeta) {
    def SAFE_ID = /[A-Za-z0-9][A-Za-z0-9._-]*/
    if (!(rawMeta instanceof Map)) {
        throw new IllegalArgumentException(
            "peak metadata must be a map, got ${rawMeta?.getClass()?.name}"
        )
    }
    def meta = new LinkedHashMap(rawMeta)
    def sampleId = meta.sample_id?.toString()
    if (sampleId == null || !sampleId.matches(SAFE_ID)) {
        throw new IllegalArgumentException(
            "peak sample_id must match ${SAFE_ID}, got ${sampleId}"
        )
    }
    def control = isControlMeta(meta)
    if (control) {
        if (meta.assay_target != null &&
            !meta.assay_target.toString().equalsIgnoreCase('IgG')) {
            throw new IllegalArgumentException(
                "control sample ${sampleId} must have assay_target IgG"
            )
        }
    } else {
        def controlId = meta.control_id?.toString()
        if (controlId == null || !controlId.matches(SAFE_ID)) {
            throw new IllegalArgumentException(
                "target ${sampleId} control_id must match ${SAFE_ID}, " +
                "got ${controlId}"
            )
        }
        meta.control_id = controlId
    }
    meta.sample_id = sampleId
    meta
}

def controlsBySampleId(controlRows) {
    def controls = new LinkedHashMap()
    controlRows.each { row ->
        def meta = row[0]
        def sampleId = meta.sample_id.toString()
        if (controls.containsKey(sampleId)) {
            throw new IllegalStateException(
                "duplicate analysis BAMs were supplied for control ${sampleId}"
            )
        }
        controls[sampleId] = row
    }
    controls
}

def validateMacsGenomeSize(rawValue) {
    def text = rawValue?.toString()
    if (text == null || text.isEmpty()) {
        throw new IllegalArgumentException(
            "macs_genome_size must be a MACS shortcut or a positive integer"
        )
    }
    if (text ==~ /[0-9]+/) {
        if (new BigInteger(text).signum() <= 0) {
            throw new IllegalArgumentException(
                "macs_genome_size integer must be positive, got ${text}"
            )
        }
        return text
    }
    if (!(text ==~ /[A-Za-z][A-Za-z0-9._-]*/)) {
        throw new IllegalArgumentException(
            "macs_genome_size shortcut contains unsafe characters: ${text}"
        )
    }
    text
}

workflow PEAKS {
    take:
    analysis_bams
    blacklist
    macs_genome_size
    motif_use_narrow_peaks

    main:
    safe_bams = analysis_bams.map { meta, bam, bai ->
        tuple(validatePeakMeta(meta), bam, bai)
    }

    controls_by_id = safe_bams
        .filter { meta, bam, bai -> isControlMeta(meta) }
        .collect(flat: false)
        .map { control_rows -> controlsBySampleId(control_rows) }

    target_bams = safe_bams.filter {
        meta, bam, bai -> !isControlMeta(meta)
    }

    paired_bams = target_bams
        .combine(controls_by_id)
        .map { meta, target_bam, target_bai, control_map ->
            def targetId = meta.sample_id.toString()
            def controlId = meta.control_id.toString()
            if (targetId == controlId) {
                throw new IllegalStateException(
                    "target ${targetId} cannot use itself as an IgG control"
                )
            }
            def control_row = control_map[controlId]
            if (control_row == null) {
                throw new IllegalStateException(
                    "no analysis BAM for control_id ${controlId} " +
                    "referenced by target ${targetId}"
                )
            }
            def control_meta = control_row[0]
            tuple(
                meta,
                target_bam,
                target_bai,
                control_meta,
                control_row[1],
                control_row[2]
            )
        }

    safe_genome_size = macs_genome_size.map { value ->
        validateMacsGenomeSize(value)
    }
    narrow_enabled = motif_use_narrow_peaks.map { enabled ->
        if (!(enabled instanceof Boolean)) {
            throw new IllegalArgumentException(
                "motif_use_narrow_peaks must be a boolean, got ${enabled}"
            )
        }
        enabled
    }

    /*
     * collect(flat: false) turns either Channel.empty() or a one-path queue
     * into a reusable value containing [] or [path]. FILTER_BLACKLIST can
     * therefore finalize every broadPeak without waiting on an empty queue.
     */
    reusable_blacklist = blacklist.collect(flat: false)
        .map { blacklist_rows ->
            if (blacklist_rows.size() > 1) {
                throw new IllegalArgumentException(
                    "PEAKS accepts at most one blacklist, got " +
                    blacklist_rows.size()
                )
            }
            blacklist_rows
        }

    MACS2_BROAD(paired_bams, safe_genome_size)
    MACS2_NARROW(paired_bams, safe_genome_size, narrow_enabled)
    FILTER_BLACKLIST(MACS2_BROAD.out.peaks, reusable_blacklist)

    logs_ch = MACS2_BROAD.out.logs.mix(
        MACS2_NARROW.out.logs,
        FILTER_BLACKLIST.out.logs
    )
    versions_ch = MACS2_BROAD.out.versions.mix(
        MACS2_NARROW.out.versions,
        FILTER_BLACKLIST.out.versions
    )

    emit:
    raw_broad_peaks = MACS2_BROAD.out.peaks
    final_broad_peaks = FILTER_BLACKLIST.out.peaks
    broad_auxiliary = MACS2_BROAD.out.auxiliary
    motif_narrow_peaks = MACS2_NARROW.out.motif_peaks
    motif_summits = MACS2_NARROW.out.motif_summits
    motif_auxiliary = MACS2_NARROW.out.auxiliary
    logs = logs_ch
    versions = versions_ch
}
