include { PREPARE_MOTIF_SEQUENCES } from '../../modules/local/prepare_motif_sequences'
include { AME } from '../../modules/local/ame'
include { STREME } from '../../modules/local/streme'
include { FIMO } from '../../modules/local/fimo'
include { MOTIF_SUMMARY } from '../../modules/local/motif_summary'

def validateMotifMeta(rawMeta, context) {
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
            "${context} is_control must be a boolean for ${sampleId}"
        )
    }
    if (meta.is_control) {
        throw new IllegalArgumentException(
            "${context} is_control must be false for motif sample ${sampleId}"
        )
    }
    def expectedMotif = meta.expected_motif?.toString()
    if (expectedMotif == null || expectedMotif.trim().isEmpty()) {
        throw new IllegalArgumentException(
            "${context} expected_motif must not be blank for ${sampleId}"
        )
    }
    try {
        java.util.regex.Pattern.compile(expectedMotif)
    } catch (java.util.regex.PatternSyntaxException error) {
        throw new IllegalArgumentException(
            "${context} expected_motif is not a valid regular expression " +
            "for ${sampleId}: ${error.description}"
        )
    }
    meta.sample_id = sampleId
    meta.expected_motif = expectedMotif
    meta
}

def uniqueMotifRows(rows, label) {
    def bySample = new LinkedHashMap()
    rows.each { row ->
        def sampleId = row[0].sample_id
        if (bySample.containsKey(sampleId)) {
            throw new IllegalStateException(
                "duplicate ${label} sample_id ${sampleId}"
            )
        }
        bySample[sampleId] = row
    }
    bySample
}

def singletonMotifPath(rows, label) {
    if (rows.size() != 1) {
        throw new IllegalArgumentException(
            "MOTIFS requires exactly one ${label}, got ${rows.size()}"
        )
    }
    rows[0]
}

def validateOptionalBlacklist(rows) {
    if (rows.size() > 1) {
        throw new IllegalArgumentException(
            "MOTIFS accepts at most one blacklist, got ${rows.size()}"
        )
    }
    rows
}

def validateNarrowSelection(rawEnabled) {
    if (!(rawEnabled instanceof Boolean)) {
        throw new IllegalArgumentException(
            "motif_use_narrow_peaks must be a boolean, got ${rawEnabled}"
        )
    }
    rawEnabled
}

def validateWorkflowMotifWindow(rawWindow) {
    def text = rawWindow?.toString()
    if (text == null || !(text ==~ /[0-9]+/)) {
        throw new IllegalArgumentException(
            "motif_window must be a positive integer, got ${text}"
        )
    }
    def window = new BigInteger(text)
    if (window.signum() <= 0 || window > Integer.MAX_VALUE) {
        throw new IllegalArgumentException(
            "motif_window must be a positive integer, got ${text}"
        )
    }
    window.toString()
}

workflow MOTIFS {
    take:
    final_motif_peaks
    motif_summits
    fasta
    blacklist
    motif_db
    motif_window
    motif_use_narrow_peaks

    main:
    safe_broad_by_sample = final_motif_peaks.map { meta, peaks ->
        tuple(validateMotifMeta(meta, 'final motif peak'), peaks)
    }.collect(flat: false).map { rows ->
        uniqueMotifRows(rows, 'final motif peak')
    }

    safe_summits_by_sample = motif_summits.map { meta, summits ->
        tuple(validateMotifMeta(meta, 'motif summit'), summits)
    }.collect(flat: false).map { rows ->
        uniqueMotifRows(rows, 'motif summit')
    }

    narrow_enabled = motif_use_narrow_peaks.map { enabled ->
        validateNarrowSelection(enabled)
    }
    safe_window = motif_window.map { window ->
        validateWorkflowMotifWindow(window)
    }
    reusable_fasta = fasta.collect(flat: false).map { rows ->
        singletonMotifPath(rows, 'reference FASTA')
    }
    reusable_motif_db = motif_db.collect(flat: false).map { rows ->
        singletonMotifPath(rows, 'motif database')
    }
    reusable_blacklist = blacklist.collect(flat: false).map { rows ->
        validateOptionalBlacklist(rows)
    }

    selected_inputs = safe_broad_by_sample
        .combine(safe_summits_by_sample)
        .combine(narrow_enabled)
        .flatMap { broadBySample, summitBySample, enabled ->
            if (enabled) {
                def missing = broadBySample.keySet() - summitBySample.keySet()
                def extra = summitBySample.keySet() - broadBySample.keySet()
                if (missing) {
                    throw new IllegalStateException(
                        "missing motif summit for sample(s): " +
                        missing.sort().join(', ')
                    )
                }
                if (extra) {
                    throw new IllegalStateException(
                        "motif summits lack final peaks for sample(s): " +
                        extra.sort().join(', ')
                    )
                }
            }
            broadBySample.values().collect { broadRow ->
                def meta = broadRow[0]
                def broadPeaks = broadRow[1]
                if (enabled) {
                    def summitRow = summitBySample[meta.sample_id]
                    if (
                        summitRow[0].sample_id != meta.sample_id ||
                        summitRow[0].expected_motif != meta.expected_motif
                    ) {
                        throw new IllegalStateException(
                            "motif summit metadata mismatch for ${meta.sample_id}"
                        )
                    }
                    tuple(meta, summitRow[1], 'summit', broadPeaks)
                } else {
                    tuple(meta, broadPeaks, 'midpoint', broadPeaks)
                }
            }
        }

    PREPARE_MOTIF_SEQUENCES(
        selected_inputs,
        reusable_fasta,
        reusable_blacklist,
        safe_window
    )
    AME(PREPARE_MOTIF_SEQUENCES.out.sequences, reusable_motif_db)
    STREME(PREPARE_MOTIF_SEQUENCES.out.sequences)
    FIMO(PREPARE_MOTIF_SEQUENCES.out.sequences, reusable_motif_db)

    keyed_ame = AME.out.results.map { meta, ameDir ->
        tuple(meta.sample_id, meta, ameDir)
    }
    keyed_fimo = FIMO.out.results.map { meta, fimoDir ->
        tuple(meta.sample_id, meta, fimoDir)
    }
    summary_inputs = keyed_ame.join(
        keyed_fimo,
        by: 0,
        failOnDuplicate: true,
        failOnMismatch: true
    ).map {
        sampleId, ameMeta, ameDir, fimoMeta, fimoDir ->
        if (
            ameMeta.sample_id != fimoMeta.sample_id ||
            ameMeta.expected_motif != fimoMeta.expected_motif
        ) {
            throw new IllegalStateException(
                "motif result metadata mismatch for ${sampleId}"
            )
        }
        tuple(ameMeta, ameDir, fimoDir)
    }

    MOTIF_SUMMARY(summary_inputs, safe_window)

    motif_metrics_ch = MOTIF_SUMMARY.out.qc.map {
        meta, json, metricsTsv, positions ->
        tuple(meta, metricsTsv)
    }
    statuses_ch = PREPARE_MOTIF_SEQUENCES.out.sequences.map {
        meta, foreground, background, status ->
        tuple(meta, 'sequence_preparation', status)
    }.mix(
        AME.out.status.map {
            meta, status -> tuple(meta, 'ame', status)
        },
        STREME.out.status.map {
            meta, status -> tuple(meta, 'streme', status)
        },
        FIMO.out.status.map {
            meta, status -> tuple(meta, 'fimo', status)
        }
    )
    logs_ch = PREPARE_MOTIF_SEQUENCES.out.logs.mix(
        AME.out.logs,
        STREME.out.logs,
        FIMO.out.logs
    )
    versions_ch = PREPARE_MOTIF_SEQUENCES.out.versions.mix(
        AME.out.versions,
        STREME.out.versions,
        FIMO.out.versions,
        MOTIF_SUMMARY.out.versions
    )

    emit:
    sequence_inputs = PREPARE_MOTIF_SEQUENCES.out.sequences
    sequence_windows = PREPARE_MOTIF_SEQUENCES.out.windows
    known_motifs = AME.out.results
    de_novo_motifs = STREME.out.results
    scans = FIMO.out.results
    expected_motif_qc = MOTIF_SUMMARY.out.qc
    motif_metrics = motif_metrics_ch
    statuses = statuses_ch
    logs = logs_ch
    versions = versions_ch
}
