nextflow.enable.dsl = 2

include { DEMULTIPLEX } from './subworkflows/local/demultiplex'
include { ALIGN_QC } from './subworkflows/local/align_qc'
include { PEAKS } from './subworkflows/local/peaks'
include { MOTIFS } from './subworkflows/local/motifs'
include { ENRICHMENT } from './subworkflows/local/enrichment'
include { QC } from './subworkflows/local/qc'

def parameterText(rawValue, label, required = true) {
    def text = rawValue == null ? null : rawValue.toString().trim()
    if (required && (text == null || text.isEmpty())) {
        throw new IllegalArgumentException("${label} must not be blank")
    }
    text == null || text.isEmpty() ? null : text
}


def resolveLocalPath(rawValue, label, launchBase) {
    def text = parameterText(rawValue, label)
    if (['*', '?', '[', ']', '{', '}'].any { marker -> text.contains(marker) }) {
        throw new IllegalArgumentException(
            "${label} must be one explicit local path, not a glob: ${text}"
        )
    }
    java.nio.file.Path path = java.nio.file.Paths.get(text)
    if (!path.isAbsolute()) {
        path = java.nio.file.Paths.get(launchBase.toString()).resolve(path)
    }
    path.normalize().toAbsolutePath()
}


def validateRegularFile(rawValue, label, launchBase, required = false) {
    def text = parameterText(rawValue, label, required)
    if (text == null) {
        return null
    }
    def path = resolveLocalPath(text, label, launchBase)
    if (
        !java.nio.file.Files.isRegularFile(path) ||
        !java.nio.file.Files.isReadable(path)
    ) {
        throw new IllegalArgumentException(
            "${label} is not a readable regular file: ${path}"
        )
    }
    path.toString()
}


def validateBowtie2IndexPrefix(rawValue, launchBase) {
    def text = parameterText(rawValue, '--bowtie2_index', false)
    if (text == null) {
        return null
    }
    def prefix = resolveLocalPath(text, '--bowtie2_index', launchBase).toString()
    def suffixes = ['1', '2', '3', '4', 'rev.1', 'rev.2']
    def small = suffixes.collect { suffix ->
        java.nio.file.Paths.get("${prefix}.${suffix}.bt2")
    }
    def large = suffixes.collect { suffix ->
        java.nio.file.Paths.get("${prefix}.${suffix}.bt2l")
    }
    def completeSmall = small.every {
        java.nio.file.Files.isRegularFile(it) &&
        java.nio.file.Files.isReadable(it)
    }
    def completeLarge = large.every {
        java.nio.file.Files.isRegularFile(it) &&
        java.nio.file.Files.isReadable(it)
    }
    if (completeSmall == completeLarge) {
        throw new IllegalArgumentException(
            "--bowtie2_index must identify exactly one complete six-file " +
            ".bt2 or .bt2l set (1,2,3,4,rev.1,rev.2): ${prefix}"
        )
    }
    prefix
}


def validateFasta(rawValue, launchBase) {
    def fasta = validateRegularFile(rawValue, '--fasta', launchBase, false)
    if (fasta == null) {
        return null
    }
    def firstNonBlank = null
    java.nio.file.Files.newBufferedReader(
        java.nio.file.Paths.get(fasta)
    ).withCloseable { reader ->
        firstNonBlank = reader.lines()
            .filter { line -> !line.trim().isEmpty() }
            .findFirst()
            .orElse(null)
    }
    if (firstNonBlank == null || !firstNonBlank.startsWith('>')) {
        throw new IllegalArgumentException(
            "--fasta must begin with a FASTA header: ${fasta}"
        )
    }
    fasta
}


def resolveManifestPath(value, manifestPathText, rowNumber, column, label) {
    def manifestPath = java.nio.file.Paths.get(manifestPathText)
    java.nio.file.Path candidate = java.nio.file.Paths.get(value)
    if (!candidate.isAbsolute()) {
        candidate = manifestPath.parent.resolve(candidate)
    }
    candidate = candidate.normalize().toAbsolutePath()
    if (
        !java.nio.file.Files.isRegularFile(candidate) ||
        !java.nio.file.Files.isReadable(candidate)
    ) {
        throw new IllegalArgumentException(
            "${label} row ${rowNumber} ${column} is not a readable regular file: ${candidate}"
        )
    }
    candidate
}


def validateChipseqManifestAtLaunch(
    rawValue,
    fastaPathText,
    launchBase,
    projectRoot
) {
    def manifestPathText = validateRegularFile(
        rawValue,
        '--chipseq_input',
        launchBase,
        false
    )
    if (manifestPathText == null) {
        return null
    }
    if (!manifestPathText.toLowerCase().endsWith('.csv')) {
        throw new IllegalArgumentException(
            "--chipseq_input must be a CSV manifest: ${manifestPathText}"
        )
    }
    if (fastaPathText == null) {
        throw new IllegalArgumentException(
            "peak enrichment requires --fasta when --chipseq_input is supplied"
        )
    }

    def validationScript = java.nio.file.Paths.get(projectRoot.toString())
        .resolve('bin/peak_enrichment.py')
        .normalize()
        .toAbsolutePath()
    if (!java.nio.file.Files.isRegularFile(validationScript)) {
        throw new IllegalArgumentException(
            "missing ChIP-seq validation utility: ${validationScript}"
        )
    }
    def command = [
        'python3',
        validationScript.toString(),
        'validate-chipseq-manifest',
        '--manifest',
        manifestPathText,
        '--fasta',
        fastaPathText,
    ]
    def validator
    try {
        validator = new ProcessBuilder(command)
            .redirectErrorStream(true)
            .start()
    } catch (java.io.IOException error) {
        throw new IllegalArgumentException(
            "python3 is required for ChIP-seq manifest launch validation",
            error
        )
    }
    def validationOutput = validator.inputStream.getText('UTF-8').trim()
    def validationStatus = validator.waitFor()
    if (validationStatus != 0) {
        throw new IllegalArgumentException(
            validationOutput ?: "ChIP-seq manifest validation failed"
        )
    }
    manifestPathText
}


def chipseqReferenceEntry(row, manifestPathText) {
    def referenceId = row.reference_id?.toString()?.trim()
    def tf = row.tf?.toString()?.trim()
    def peakFile = row.peak_file?.toString()?.trim()
    if (referenceId == null || tf == null || peakFile == null) {
        throw new IllegalArgumentException(
            "--chipseq_input rows require reference_id, tf, and peak_file"
        )
    }
    def peakPath = resolveManifestPath(
        peakFile,
        manifestPathText,
        'parsed',
        'peak_file',
        '--chipseq_input'
    )
    tuple(
        [
            reference_id: referenceId,
            tf: tf,
            reference_type: 'chipseq',
        ],
        peakPath
    )
}


def validateMemeDatabase(rawValue, launchBase) {
    def motifDb = validateRegularFile(
        rawValue,
        '--motif_db',
        launchBase,
        false
    )
    if (motifDb == null) {
        return null
    }
    def firstNonBlank = null
    def hasMotif = false
    java.nio.file.Files.newBufferedReader(
        java.nio.file.Paths.get(motifDb)
    ).withCloseable { reader ->
        reader.eachLine { line ->
            if (firstNonBlank == null && !line.trim().isEmpty()) {
                firstNonBlank = line
            }
            if (line ==~ /\s*MOTIF\s+\S+(?:\s+\S+)?\s*/) {
                hasMotif = true
            }
        }
    }
    if (
        firstNonBlank == null ||
        !(firstNonBlank ==~ /MEME version\s+[0-9]+(?:\.[0-9]+)*/)
    ) {
        throw new IllegalArgumentException(
            "--motif_db must be a MEME-format database beginning with " +
            "'MEME version': ${motifDb}"
        )
    }
    if (!hasMotif) {
        throw new IllegalArgumentException(
            "--motif_db contains no MOTIF records: ${motifDb}"
        )
    }
    motifDb
}


def validateNonNegativeInteger(rawValue, label) {
    def text = parameterText(rawValue, label)
    if (!(text ==~ /[0-9]+/)) {
        throw new IllegalArgumentException(
            "${label} must be a non-negative integer, got ${text}"
        )
    }
    def value = new BigInteger(text)
    if (value > Integer.MAX_VALUE) {
        throw new IllegalArgumentException(
            "${label} is outside the supported integer range: ${text}"
        )
    }
    value.intValue()
}


def validatePositiveInteger(rawValue, label) {
    def value = validateNonNegativeInteger(rawValue, label)
    if (value == 0) {
        throw new IllegalArgumentException(
            "${label} must be a positive integer, got 0"
        )
    }
    value
}


def validateInteger(rawValue, label) {
    def text = parameterText(rawValue, label)
    if (!(text ==~ /-?[0-9]+/)) {
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
    value.intValue()
}


def validateUnitInterval(rawValue, label) {
    def text = parameterText(rawValue, label)
    def value
    try {
        value = new BigDecimal(text)
    } catch (NumberFormatException error) {
        throw new IllegalArgumentException(
            "${label} must be a number in the interval (0, 1], got ${text}"
        )
    }
    if (value <= BigDecimal.ZERO || value > BigDecimal.ONE) {
        throw new IllegalArgumentException(
            "${label} must be a number in the interval (0, 1], got ${text}"
        )
    }
    value.doubleValue()
}


def validateMapq(rawValue) {
    def value = validateNonNegativeInteger(rawValue, '--min_mapq')
    if (value > 255) {
        throw new IllegalArgumentException(
            "--min_mapq must be between 0 and 255, got ${value}"
        )
    }
    value
}


def validateBoolean(rawValue, label) {
    if (!(rawValue instanceof Boolean)) {
        throw new IllegalArgumentException(
            "${label} must be true or false, got ${rawValue}"
        )
    }
    rawValue
}


def validateMacsGenomeSize(rawValue) {
    def text = parameterText(rawValue, '--macs_genome_size')
    if (text ==~ /[0-9]+/) {
        if (new BigInteger(text).signum() <= 0) {
            throw new IllegalArgumentException(
                "--macs_genome_size integer must be positive"
            )
        }
        return text
    }
    if (!(text ==~ /[A-Za-z][A-Za-z0-9._-]*/)) {
        throw new IllegalArgumentException(
            "--macs_genome_size must be a positive integer or path-safe " +
            "MACS shortcut, got ${text}"
        )
    }
    text
}


def validateFixedSetting(rawValue, expectedValue, label) {
    if (rawValue.toString() != expectedValue.toString()) {
        throw new IllegalArgumentException(
            "${label} must remain ${expectedValue} for NanoScope compatibility, " +
            "got ${rawValue}"
        )
    }
    expectedValue
}


def validateOutputDirectory(rawValue, launchBase) {
    def path = resolveLocalPath(rawValue, '--outdir', launchBase)
    if (
        java.nio.file.Files.exists(path) &&
        !java.nio.file.Files.isDirectory(path)
    ) {
        throw new IllegalArgumentException(
            "--outdir exists but is not a directory: ${path}"
        )
    }
    path.toString()
}


def validatePipelineParameters(rawParams, launchBase, projectRoot) {
    def validated = new LinkedHashMap()
    validated.input = validateRegularFile(
        rawParams.input,
        '--input',
        launchBase,
        true
    )
    if (!validated.input.toLowerCase().endsWith('.csv')) {
        throw new IllegalArgumentException(
            "--input must be a CSV manifest: ${validated.input}"
        )
    }
    validated.fasta = validateFasta(rawParams.fasta, launchBase)
    validated.bowtie2_index = validateBowtie2IndexPrefix(
        rawParams.bowtie2_index,
        launchBase
    )
    if (validated.fasta == null && validated.bowtie2_index == null) {
        throw new IllegalArgumentException(
            "pipeline requires --fasta or --bowtie2_index"
        )
    }
    validated.blacklist = validateRegularFile(
        rawParams.blacklist,
        '--blacklist',
        launchBase,
        false
    )
    validated.gtf = validateRegularFile(
        rawParams.gtf,
        '--gtf',
        launchBase,
        false
    )
    validated.tss_bed = validateRegularFile(
        rawParams.tss_bed,
        '--tss_bed',
        launchBase,
        false
    )
    validated.motif_db = validateMemeDatabase(
        rawParams.motif_db,
        launchBase
    )
    if (validated.motif_db != null && validated.fasta == null) {
        throw new IllegalArgumentException(
            "motif analysis requires --fasta when --motif_db is supplied"
        )
    }
    validated.chipseq_input = validateChipseqManifestAtLaunch(
        rawParams.chipseq_input,
        validated.fasta,
        launchBase,
        projectRoot
    )
    validated.outdir = validateOutputDirectory(rawParams.outdir, launchBase)
    validated.barcode_mismatches = validateNonNegativeInteger(
        rawParams.barcode_mismatches,
        '--barcode_mismatches'
    )
    validated.allow_empty = validateBoolean(
        rawParams.allow_empty,
        '--allow_empty'
    )
    validated.min_mapq = validateMapq(rawParams.min_mapq)
    validated.macs_genome_size = validateMacsGenomeSize(
        rawParams.macs_genome_size
    )
    validated.macs_llocal = validateFixedSetting(
        rawParams.macs_llocal,
        100000,
        '--macs_llocal'
    )
    validated.macs_keep_dup = validateFixedSetting(
        rawParams.macs_keep_dup,
        1,
        '--macs_keep_dup'
    )
    validated.macs_broad_cutoff = validateFixedSetting(
        rawParams.macs_broad_cutoff,
        0.1,
        '--macs_broad_cutoff'
    )
    validated.macs_max_gap = validateFixedSetting(
        rawParams.macs_max_gap,
        1000,
        '--macs_max_gap'
    )
    validated.motif_use_narrow_peaks = validateBoolean(
        rawParams.motif_use_narrow_peaks,
        '--motif_use_narrow_peaks'
    )
    validated.motif_window = validatePositiveInteger(
        rawParams.motif_window,
        '--motif_window'
    )
    validated.enrichment_permutations = rawParams.enrichment_permutations == null
        ? 1000
        : validatePositiveInteger(
            rawParams.enrichment_permutations,
            '--enrichment_permutations'
        )
    validated.enrichment_seed = rawParams.enrichment_seed == null
        ? 1729
        : validateInteger(
            rawParams.enrichment_seed,
            '--enrichment_seed'
        )
    validated.enrichment_gc_tolerance =
        rawParams.enrichment_gc_tolerance == null
            ? 0.02d
            : validateUnitInterval(
                rawParams.enrichment_gc_tolerance,
                '--enrichment_gc_tolerance'
            )
    validated
}


def optionalPathChannel(pathText) {
    pathText == null
        ? Channel.empty()
        : Channel.fromPath(pathText, checkIfExists: true)
}


def versionPath(versionRecord) {
    if (!(versionRecord instanceof java.util.Collection)) {
        return versionRecord
    }
    if (versionRecord.isEmpty()) {
        throw new IllegalArgumentException(
            "version output contains an empty tuple"
        )
    }
    versionPath(versionRecord.last())
}


process WRITE_PIPELINE_PARAMETERS {
    tag 'validated parameters'
    label 'process_light'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy',
        overwrite: true

    input:
    val encoded_parameters

    output:
    path "validated_parameters.json", emit: parameters
    path "pipeline_parameters_versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    printf '%s' "${encoded_parameters}" | python -c \
        'import base64, sys; sys.stdout.buffer.write(base64.b64decode(sys.stdin.buffer.read()))' \
        > "validated_parameters.json"
    printf 'WRITE_PIPELINE_PARAMETERS:\\n  python: ' \
        > "pipeline_parameters_versions.yml"
    python --version 2>&1 >> "pipeline_parameters_versions.yml"
    """
}


process COLLECT_VERSIONS {
    tag 'software versions'
    label 'process_light'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy',
        overwrite: true

    input:
    path version_files, stageAs: 'versions??/*', arity: '1..*'

    output:
    path "software_versions.yml", emit: versions

    script:
    '''
    set -euo pipefail
    python <<'PY'
from pathlib import Path

files = sorted(Path(".").glob("versions*/*"))
if not files:
    raise SystemExit("no version records were supplied")
with open("software_versions.yml", "w", encoding="utf-8") as output:
    output.write("pipeline_software_versions:\\n")
    for index, path in enumerate(files, start=1):
        output.write(f"  record_{index:04d}:\\n")
        output.write(f"    source: {path.name}\\n")
        output.write("    content: |\\n")
        for line in path.read_text(encoding="utf-8").splitlines():
            output.write(f"      {line}\\n")
    output.write("  collector:\\n")
    output.write("    implementation: repository\\n")
PY
    '''
}


process WRITE_COMPLETION_SUMMARY {
    tag 'completion summary'
    label 'process_light'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy',
        overwrite: true

    input:
    path multiqc_report, stageAs: 'report/multiqc_report.html'
    path combined_summary, stageAs: 'report/combined_target_qc.tsv'
    path enrichment_files, stageAs: 'enrichment??/*'
    path software_versions, stageAs: 'pipeline/software_versions.yml'
    path validated_parameters, stageAs: 'pipeline/validated_parameters.json'

    output:
    path "run_summary.txt", emit: summary

    script:
    '''
    set -euo pipefail
    python <<'PY'
import csv
import json
from pathlib import Path

with open("pipeline/validated_parameters.json", encoding="utf-8") as handle:
    parameters = json.load(handle)
with open(
    "report/combined_target_qc.tsv",
    encoding="utf-8",
    newline="",
) as handle:
    targets = list(csv.DictReader(handle, delimiter="\\t"))
if not Path("report/multiqc_report.html").stat().st_size:
    raise SystemExit("MultiQC report is empty")
if not Path("pipeline/software_versions.yml").stat().st_size:
    raise SystemExit("software version manifest is empty")
enrichment_tables = sorted(Path().glob("enrichment*/*"))
peak_enrichment_tables = [
    path for path in enrichment_tables if path.name == "peak_enrichment.tsv"
]
if len(peak_enrichment_tables) > 1:
    raise SystemExit("expected at most one peak_enrichment.tsv in completion summary inputs")
enrichment_by_name = {}
for path in enrichment_tables:
    if path.name in enrichment_by_name:
        raise SystemExit(f"duplicate enrichment summary input: {path.name}")
    enrichment_by_name[path.name] = path
with open("run_summary.txt", "w", encoding="utf-8") as output:
    output.write("status\\tsucceeded\\n")
    output.write(f"target_count\\t{len(targets)}\\n")
    output.write(f"motif_enabled\\t{str(bool(parameters['motif_db'])).lower()}\\n")
    output.write("report\\treports/multiqc/multiqc_report.html\\n")
    output.write("target_summary\\treports/summary/combined_target_qc.tsv\\n")
    if peak_enrichment_tables:
        required_plots = [
            "observed_vs_null.png",
            "matrix_random.png",
            "matrix_length_matched.png",
            "matrix_gc_matched.png",
            "matrix_length_gc_matched.png",
        ]
        missing_plots = [
            name for name in required_plots if name not in enrichment_by_name
        ]
        if missing_plots:
            raise SystemExit(
                "missing enrichment plots for completion summary: "
                + ", ".join(missing_plots)
            )
        output.write("peak_enrichment\\tenrichment/peak_enrichment.tsv\\n")
        output.write("enrichment_plot\\tenrichment/plots/observed_vs_null.png\\n")
        for model in (
            "random",
            "length_matched",
            "gc_matched",
            "length_gc_matched",
        ):
            output.write(
                f"enrichment_heatmap_{model}\\t"
                f"enrichment/plots/matrix_{model}.png\\n"
            )
    output.write("versions\\tpipeline_info/software_versions.yml\\n")
PY
    '''
}


workflow NANOCUT {
    main:
    validated = validatePipelineParameters(params, launchDir, projectDir)
    params.outdir = validated.outdir

    manifest_ch = Channel.fromPath(validated.input, checkIfExists: true)
    blacklist_ch = optionalPathChannel(validated.blacklist)
    gtf_ch = optionalPathChannel(validated.gtf)
    tss_bed_ch = optionalPathChannel(validated.tss_bed)

    barcode_mismatches_ch = Channel.value(validated.barcode_mismatches)
    allow_empty_ch = Channel.value(validated.allow_empty)
    min_mapq_ch = Channel.value(validated.min_mapq)
    macs_genome_size_ch = Channel.value(validated.macs_genome_size)
    narrow_peaks_ch = Channel.value(
        validated.motif_db != null && validated.motif_use_narrow_peaks
    )
    enrichment_results_ch = Channel.empty()
    enrichment_status_ch = Channel.empty()
    enrichment_plots_ch = Channel.empty()
    enrichment_dashboard_files_ch = Channel.empty()
    enrichment_versions_ch = Channel.empty()
    chipseq_reference_inputs_ch = Channel.empty()

    DEMULTIPLEX(
        manifest_ch,
        barcode_mismatches_ch,
        allow_empty_ch
    )
    ALIGN_QC(
        DEMULTIPLEX.out.reads,
        validated.fasta,
        validated.bowtie2_index,
        min_mapq_ch
    )
    PEAKS(
        ALIGN_QC.out.analysis_bam,
        blacklist_ch,
        macs_genome_size_ch,
        narrow_peaks_ch
    )

    motif_metrics_ch = Channel.empty()
    motif_versions_ch = Channel.empty()
    if (validated.motif_db != null) {
        motif_fasta_ch = Channel.fromPath(
            validated.fasta,
            checkIfExists: true
        )
        motif_db_ch = Channel.fromPath(
            validated.motif_db,
            checkIfExists: true
        )
        motif_blacklist_ch = optionalPathChannel(validated.blacklist)
        motif_use_narrow_peaks_ch = Channel.value(
            validated.motif_use_narrow_peaks
        )
        motif_window_ch = Channel.value(validated.motif_window)
        MOTIFS(
            PEAKS.out.final_broad_peaks,
            PEAKS.out.motif_summits,
            motif_fasta_ch,
            motif_blacklist_ch,
            motif_db_ch,
            motif_window_ch,
            motif_use_narrow_peaks_ch
        )
        motif_metrics_ch = MOTIFS.out.motif_metrics
        motif_versions_ch = MOTIFS.out.versions
    }

    if (validated.chipseq_input != null) {
        chipseq_reference_inputs_ch = Channel
            .fromPath(validated.chipseq_input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                chipseqReferenceEntry(row, validated.chipseq_input)
            }
            .collect(flat: false)
            .map { entries ->
                if (entries.isEmpty()) {
                    throw new IllegalArgumentException(
                        "--chipseq_input manifest contains no rows"
                    )
                }
                def orderedEntries = entries.sort { left, right ->
                    left[0].reference_id <=> right[0].reference_id
                }
                def referenceRows = orderedEntries.collect { entry ->
                    entry[0]
                }
                def encodedRows = groovy.json.JsonOutput
                    .toJson(referenceRows)
                    .bytes
                    .encodeBase64()
                    .toString()
                tuple(
                    encodedRows,
                    orderedEntries.collect { entry -> entry[1] }
                )
            }
        enrichment_fasta_ch = Channel.fromPath(
            validated.fasta,
            checkIfExists: true
        )
        enrichment_blacklist_ch = optionalPathChannel(validated.blacklist)
        enrichment_permutations_ch = Channel.value(
            validated.enrichment_permutations
        )
        enrichment_seed_ch = Channel.value(validated.enrichment_seed)
        enrichment_gc_tolerance_ch = Channel.value(
            validated.enrichment_gc_tolerance
        )
        analysis_metadata_ch = PEAKS.out.final_broad_peaks.map {
            meta, peaks -> meta
        }

        ENRICHMENT(
            PEAKS.out.final_broad_peaks,
            analysis_metadata_ch,
            chipseq_reference_inputs_ch,
            enrichment_fasta_ch,
            enrichment_blacklist_ch,
            enrichment_permutations_ch,
            enrichment_seed_ch,
            enrichment_gc_tolerance_ch
        )
        enrichment_results_ch = ENRICHMENT.out.results
        enrichment_status_ch = ENRICHMENT.out.status
        enrichment_plots_ch = ENRICHMENT.out.plots
        enrichment_dashboard_files_ch = ENRICHMENT.out.results.mix(
            ENRICHMENT.out.plots
        )
        enrichment_versions_ch = ENRICHMENT.out.versions
    }

    QC(
        ALIGN_QC.out.filtered_bam,
        PEAKS.out.final_broad_peaks,
        ALIGN_QC.out.coverage,
        ALIGN_QC.out.metrics,
        DEMULTIPLEX.out.metrics,
        DEMULTIPLEX.out.fastqc,
        motif_metrics_ch,
        enrichment_dashboard_files_ch,
        gtf_ch,
        tss_bed_ch
    )

    parameter_map = new LinkedHashMap(validated)
    parameter_map.profile = workflow.profile
    parameter_json = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson(parameter_map)
    ) + '\n'
    encoded_parameters = parameter_json.bytes.encodeBase64().toString()
    WRITE_PIPELINE_PARAMETERS(Channel.value(encoded_parameters))

    all_versions_ch = DEMULTIPLEX.out.versions.mix(
        ALIGN_QC.out.versions,
        PEAKS.out.versions,
        motif_versions_ch,
        enrichment_versions_ch,
        QC.out.versions,
        WRITE_PIPELINE_PARAMETERS.out.versions
    )
    version_files_ch = all_versions_ch.map { record ->
        versionPath(record)
    }.collect()
    COLLECT_VERSIONS(version_files_ch)
    completion_summary_enrichment = QC.out.enrichment_table
        .mix(QC.out.enrichment_plot, QC.out.enrichment_heatmaps)
        .collect()
        .ifEmpty { ignored -> [] }

    WRITE_COMPLETION_SUMMARY(
        QC.out.multiqc_report,
        QC.out.combined_summary,
        completion_summary_enrichment,
        COLLECT_VERSIONS.out.versions,
        WRITE_PIPELINE_PARAMETERS.out.parameters
    )

    all_metrics_ch = DEMULTIPLEX.out.metrics.mix(
        ALIGN_QC.out.metrics,
        QC.out.target_qc,
        motif_metrics_ch
    )

    emit:
    demultiplexed_reads = DEMULTIPLEX.out.reads
    demultiplex_metrics = DEMULTIPLEX.out.metrics
    fastqc = DEMULTIPLEX.out.fastqc
    analysis_bams = ALIGN_QC.out.analysis_bam
    filtered_bams = ALIGN_QC.out.filtered_bam
    library_metrics = ALIGN_QC.out.metrics
    coverage = ALIGN_QC.out.coverage
    final_broad_peaks = PEAKS.out.final_broad_peaks
    motif_metrics = motif_metrics_ch
    enrichment_results = enrichment_results_ch
    enrichment_status = enrichment_status_ch
    enrichment_plots = enrichment_plots_ch
    target_qc = QC.out.target_qc
    combined_summary = QC.out.combined_summary
    multiqc_report = QC.out.multiqc_report
    metrics = all_metrics_ch
    versions = COLLECT_VERSIONS.out.versions
    completion_summary = WRITE_COMPLETION_SUMMARY.out.summary
}


workflow {
    NANOCUT()
}
