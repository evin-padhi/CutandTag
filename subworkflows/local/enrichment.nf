include { VALIDATE_CHIPSEQ_MANIFEST } from '../../modules/local/validate_chipseq_manifest'

def ENRICHMENT_SAFE_ID = /[A-Za-z0-9][A-Za-z0-9._-]*/

def validateEnrichmentMeta(rawMeta, context) {
    if (!(rawMeta instanceof Map)) {
        throw new IllegalArgumentException(
            "${context} metadata must be a map"
        )
    }
    def meta = new LinkedHashMap(rawMeta)
    def sampleId = meta.sample_id?.toString()
    if (sampleId == null || !sampleId.matches(ENRICHMENT_SAFE_ID)) {
        throw new IllegalArgumentException(
            "${context} sample_id must match ${ENRICHMENT_SAFE_ID}, got ${sampleId}"
        )
    }
    if (!(meta.is_control instanceof Boolean)) {
        throw new IllegalArgumentException(
            "${context} is_control must be a boolean for ${sampleId}"
        )
    }
    def assayTarget = meta.assay_target?.toString()?.trim()
    if (assayTarget == null || assayTarget.isEmpty()) {
        throw new IllegalArgumentException(
            "${context} assay_target must not be blank for ${sampleId}"
        )
    }
    meta.sample_id = sampleId
    meta.assay_target = assayTarget
    meta
}

def analysisBySampleId(rows) {
    def bySampleId = new LinkedHashMap()
    rows.each { meta ->
        def sampleId = meta.sample_id.toString()
        if (bySampleId.containsKey(sampleId)) {
            throw new IllegalStateException(
                "duplicate analysis metadata for sample ${sampleId}"
            )
        }
        bySampleId[sampleId] = meta
    }
    bySampleId
}

def enrichmentDescriptor(rawTf) {
    def tf = rawTf?.toString()?.trim()
    if (tf == null || tf.isEmpty()) {
        throw new IllegalArgumentException(
            "foreground TF must not be blank"
        )
    }
    def token = tf.replaceAll(/[^A-Za-z0-9._-]+/, '_')
    if (!token || !(token[0] ==~ /[A-Za-z0-9]/)) {
        token = "tf_${token}"
    }
    def foregroundId = "called_tf_${token}"
    if (!foregroundId.matches(ENRICHMENT_SAFE_ID)) {
        throw new IllegalArgumentException(
            "foreground_id must match ${ENRICHMENT_SAFE_ID}, got ${foregroundId}"
        )
    }
    [
        foreground_id: foregroundId,
        foreground_tf: tf,
    ]
}

def singletonEnrichmentPath(rows, label) {
    if (rows.size() != 1) {
        throw new IllegalArgumentException(
            "ENRICHMENT requires exactly one ${label}, got ${rows.size()}"
        )
    }
    rows[0]
}

def validateOptionalBlacklist(rows) {
    if (rows.size() > 1) {
        throw new IllegalArgumentException(
            "ENRICHMENT accepts at most one blacklist, got ${rows.size()}"
        )
    }
    rows
}

def validateWorkflowPositiveInteger(rawValue, label) {
    def text = rawValue?.toString()
    if (text == null || !(text ==~ /[0-9]+/)) {
        throw new IllegalArgumentException(
            "${label} must be a positive integer, got ${text}"
        )
    }
    def value = new BigInteger(text)
    if (value.signum() <= 0 || value > Integer.MAX_VALUE) {
        throw new IllegalArgumentException(
            "${label} must be a positive integer, got ${text}"
        )
    }
    value.toString()
}

def validateWorkflowInteger(rawValue, label) {
    def text = rawValue?.toString()
    if (text == null || !(text ==~ /-?[0-9]+/)) {
        throw new IllegalArgumentException(
            "${label} must be an integer, got ${text}"
        )
    }
    def value = new BigInteger(text)
    if (value < Integer.MIN_VALUE || value > Integer.MAX_VALUE) {
        throw new IllegalArgumentException(
            "${label} is outside the supported integer range: ${text}"
        )
    }
    value.toString()
}

def validateWorkflowGcTolerance(rawValue) {
    def text = rawValue?.toString()
    if (text == null) {
        throw new IllegalArgumentException(
            "gc_tolerance must be a number in the interval (0, 1]"
        )
    }
    def value
    try {
        value = new BigDecimal(text)
    } catch (NumberFormatException error) {
        throw new IllegalArgumentException(
            "gc_tolerance must be a number in the interval (0, 1], got ${text}"
        )
    }
    if (value <= BigDecimal.ZERO || value > BigDecimal.ONE) {
        throw new IllegalArgumentException(
            "gc_tolerance must be a number in the interval (0, 1], got ${text}"
        )
    }
    value.toPlainString()
}

process MERGE_FOREGROUND_PEAKS {
    tag "${foreground_meta.foreground_id}"
    label 'process_light'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    input:
    tuple val(foreground_meta),
        path(peak_files, stageAs: 'foreground??/*.broadPeak', arity: '1..*')

    output:
    tuple val(foreground_meta),
        path("${foreground_meta.foreground_id}.bed"),
        emit: foregrounds
    tuple val(foreground_meta),
        path("merge_foreground_peaks_versions.yml"),
        emit: versions

    script:
    def foregroundId = foreground_meta.foreground_id?.toString()
    def foregroundTf = foreground_meta.foreground_tf?.toString()
    if (foregroundId == null || !foregroundId.matches(ENRICHMENT_SAFE_ID)) {
        throw new IllegalArgumentException(
            "foreground_id must match ${ENRICHMENT_SAFE_ID}, got ${foregroundId}"
        )
    }
    if (foregroundTf == null || foregroundTf.trim().isEmpty()) {
        throw new IllegalArgumentException(
            "foreground_tf must not be blank for ${foregroundId}"
        )
    }

    """
    set -euo pipefail
    awk -F '\\t' '
        BEGIN { OFS = "\\t" }
        /^#/ || NF < 3 { next }
        \$2 ~ /^[0-9]+\$/ && \$3 ~ /^[0-9]+\$/ && \$3 > \$2 {
            print \$1, \$2, \$3
        }
    ' foreground*/*.broadPeak \
        | LC_ALL=C sort -t \$'\\t' -k1,1 -k2,2n -k3,3n -u \
        > "${foregroundId}.bed"

    if [[ ! -f "${foregroundId}.bed" ]]; then
        : > "${foregroundId}.bed"
    fi

    printf 'MERGE_FOREGROUND_PEAKS:\\n  implementation: repository\\n' \
        > "merge_foreground_peaks_versions.yml"
    """
}

process RUN_PEAK_ENRICHMENT {
    tag 'peak-enrichment'
    label 'process_heavy'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/peak_enrichment",
        mode: 'copy',
        overwrite: true

    input:
    val foreground_rows_b64
    path foreground_peak_files, stageAs: 'foregrounds/*', arity: '1..*'
    path external_manifest, stageAs: 'chipseq/validated.tsv'
    path fasta, stageAs: 'reference.fa'
    path blacklist_files, stageAs: 'blacklist/regions*.bed'
    val permutations
    val seed
    val gc_tolerance

    output:
    path "peak_enrichment.tsv", emit: results
    path "enrichment_status.tsv", emit: status
    path "plots", emit: plots
    path "peak_enrichment_versions.yml", emit: versions

    script:
    def safePermutations = validateWorkflowPositiveInteger(
        permutations,
        'enrichment permutations'
    )
    def safeSeed = validateWorkflowInteger(seed, 'enrichment seed')
    def safeGcTolerance = validateWorkflowGcTolerance(gc_tolerance)
    if (blacklist_files.size() > 1) {
        throw new IllegalArgumentException(
            "RUN_PEAK_ENRICHMENT accepts at most one blacklist, got ${blacklist_files.size()}"
        )
    }
    def blacklistArg = blacklist_files
        ? '--blacklist "blacklist/regions.bed"'
        : ''
    def projectBin = "${projectDir}/bin".toString()

    """
    set -euo pipefail
    printf '%s' "${foreground_rows_b64}" | python -c \
        'import base64, sys; sys.stdout.buffer.write(base64.b64decode(sys.stdin.buffer.read()))' \
        > "foreground_rows.json"

    python <<'PY'
import csv
import json
from pathlib import Path

rows = json.loads(Path("foreground_rows.json").read_text(encoding="utf-8"))

with open("foreground_manifest.tsv", "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["foreground_id", "foreground_tf", "peak_file"],
        delimiter="\\t",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

with open("chipseq/validated.tsv", encoding="utf-8", newline="") as handle:
    external_rows = list(csv.DictReader(handle, delimiter="\\t"))

with open("reference_manifest.tsv", "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["reference_id", "tf", "reference_type", "peak_file"],
        delimiter="\\t",
    )
    writer.writeheader()
    for row in external_rows:
        writer.writerow(row)
    for row in rows:
        writer.writerow(
            {
                "reference_id": row["foreground_id"],
                "tf": row["foreground_tf"],
                "reference_type": "called_tf",
                "peak_file": row["peak_file"],
            }
        )
PY

    python "${projectBin}/peak_enrichment.py" \
        --foreground-manifest "foreground_manifest.tsv" \
        --reference-manifest "reference_manifest.tsv" \
        --fasta "reference.fa" \
        --outdir "." \
        --permutations "${safePermutations}" \
        --seed "${safeSeed}" \
        --gc-tolerance "${safeGcTolerance}" \
        ${blacklistArg}

    mkdir -p "plots"
    mv "observed_vs_null.png" "plots/"
    mv matrix_*.tsv matrix_*.png "plots/"

    printf 'RUN_PEAK_ENRICHMENT:\\n  python: ' \
        > "peak_enrichment_versions.yml"
    python --version 2>&1 >> "peak_enrichment_versions.yml"
    printf '  peak_enrichment.py: repository\\n' \
        >> "peak_enrichment_versions.yml"
    """
}

workflow ENRICHMENT {
    take:
    final_broad_peaks
    analysis_metadata
    chipseq_manifest
    fasta
    blacklist
    permutations
    seed
    gc_tolerance

    main:
    safe_analysis_by_sample = analysis_metadata
        .map { meta -> validateEnrichmentMeta(meta, 'analysis metadata') }
        .collect(flat: false)
        .map { rows -> analysisBySampleId(rows) }

    reusable_chipseq_manifest = chipseq_manifest
        .collect(flat: false)
        .map { rows -> singletonEnrichmentPath(rows, 'ChIP-seq manifest') }
    reusable_fasta = fasta
        .collect(flat: false)
        .map { rows -> singletonEnrichmentPath(rows, 'reference FASTA') }
    reusable_blacklist = blacklist
        .reduce([files: []]) { holder, blacklist_path ->
            [files: holder.files + [blacklist_path]]
        }
        .map { holder -> validateOptionalBlacklist(holder.files) }
    safe_permutations = permutations.map { value ->
        validateWorkflowPositiveInteger(value, 'enrichment permutations')
    }
    safe_seed = seed.map { value ->
        validateWorkflowInteger(value, 'enrichment seed')
    }
    safe_gc_tolerance = gc_tolerance.map { value ->
        validateWorkflowGcTolerance(value)
    }

    VALIDATE_CHIPSEQ_MANIFEST(reusable_chipseq_manifest, reusable_fasta)

    grouped_foregrounds = final_broad_peaks
        .map { meta, peaks ->
            tuple(validateEnrichmentMeta(meta, 'final broad peak'), peaks)
        }
        .filter { meta, peaks -> !meta.is_control }
        .combine(safe_analysis_by_sample)
        .map { peak_meta, peak_file, analysis_by_sample ->
            def analysis_meta = analysis_by_sample[peak_meta.sample_id]
            if (analysis_meta == null) {
                throw new IllegalStateException(
                    "final broad peaks were emitted for unknown sample ${peak_meta.sample_id}"
                )
            }
            if (analysis_meta.assay_target != peak_meta.assay_target) {
                throw new IllegalStateException(
                    "assay_target mismatch for sample ${peak_meta.sample_id}"
                )
            }
            tuple(enrichmentDescriptor(analysis_meta.assay_target), peak_file)
        }
        .groupTuple(by: 0)

    MERGE_FOREGROUND_PEAKS(grouped_foregrounds)

    MERGE_FOREGROUND_PEAKS.out.foregrounds.into {
        merged_foregrounds_for_rows
        merged_foregrounds_for_paths
    }

    foreground_rows_b64 = merged_foregrounds_for_rows
        .map { foreground_meta, peak_file ->
            [
                foreground_id: foreground_meta.foreground_id,
                foreground_tf: foreground_meta.foreground_tf,
                peak_file: "foregrounds/${peak_file.name}",
            ]
        }
        .collect(flat: false)
        .map { rows ->
            groovy.json.JsonOutput.toJson(rows).bytes
                .encodeBase64()
                .toString()
        }
    foreground_peak_files = merged_foregrounds_for_paths
        .map { foreground_meta, peak_file -> peak_file }
        .collect(flat: false)

    RUN_PEAK_ENRICHMENT(
        foreground_rows_b64,
        foreground_peak_files,
        VALIDATE_CHIPSEQ_MANIFEST.out.normalized,
        reusable_fasta,
        reusable_blacklist,
        safe_permutations,
        safe_seed,
        safe_gc_tolerance
    )

    versions_ch = VALIDATE_CHIPSEQ_MANIFEST.out.versions.mix(
        MERGE_FOREGROUND_PEAKS.out.versions,
        RUN_PEAK_ENRICHMENT.out.versions
    )

    emit:
    results = RUN_PEAK_ENRICHMENT.out.results
    plots = RUN_PEAK_ENRICHMENT.out.plots
    status = RUN_PEAK_ENRICHMENT.out.status
    versions = versions_ch
}
