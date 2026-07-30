def tssSampleId(meta) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    def sampleId = meta instanceof Map ? meta.sample_id?.toString() : null
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "TSS sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    sampleId
}

process TSS_ENRICHMENT {
    tag "${meta.sample_id}"
    label 'process_heavy'

    conda "${projectDir}/envs/deeptools.yml"
    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'

    publishDir "${params.outdir}/qc/tss/${meta.sample_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(bigwig, stageAs: 'coverage.RPKM.bw'),
        val(annotation_mode),
        path(annotation_files, stageAs: 'annotation/source.annotation', arity: '1')

    output:
    tuple val(meta),
        path("*.tss.bed"),
        path("*.tss_matrix.gz"),
        path("*.tss_matrix.tsv"),
        path("*.tss_profile.png"),
        path("*.tss_profile.tsv"),
        path("*.tss_status.tsv"),
        emit: profiles
    tuple val(meta), path("tss_enrichment_versions.yml"), emit: versions

    script:
    def sampleId = tssSampleId(meta)
    if (!(annotation_mode in ['gtf', 'bed'])) {
        throw new IllegalArgumentException(
            "TSS annotation mode must be gtf or bed, got ${annotation_mode}"
        )
    }
    def outputStem = sampleId
    def prepareTss = annotation_mode == 'gtf' ? '''
    awk -F '\t' '
        BEGIN { OFS = "\t" }
        $0 !~ /^#/ && $3 == "transcript" && ($7 == "+" || $7 == "-") {
            if ($7 == "+") {
                start = $4 - 1
            } else if ($7 == "-") {
                start = $5 - 1
            }
            if (start < 0) {
                start = 0
            }
            print $1, start, start + 1, "tss_" NR, 0, $7
        }
    ' "annotation/source.annotation" \
        | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k6,6 \
        > "__OUTPUT_STEM__.tss.bed"
    if [[ ! -s "__OUTPUT_STEM__.tss.bed" ]]; then
        printf '%s\n' \
            'GTF contains no strand-bearing transcript features for TSS enrichment' \
            >&2
        exit 2
    fi
    ''' : '''
    awk -F '\t' '
        BEGIN { OFS = "\t" }
        /^#/ || NF == 0 { next }
        NF < 6 {
            print "TSS BED must be BED6 with strand in column 6" > "/dev/stderr"
            exit 2
        }
        $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= $2 {
            print "TSS BED coordinates must be non-negative half-open intervals" \
                > "/dev/stderr"
            exit 2
        }
        $6 != "+" && $6 != "-" {
            print "TSS BED6 strand must be + or -" > "/dev/stderr"
            exit 2
        }
        { print $1, $2, $3, $4, $5, $6 }
    ' "annotation/source.annotation" \
        | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k6,6 \
        > "__OUTPUT_STEM__.tss.bed"
    if [[ ! -s "__OUTPUT_STEM__.tss.bed" ]]; then
        printf '%s\n' 'TSS BED contains no usable BED6 records' >&2
        exit 2
    fi
    '''
    prepareTss = prepareTss.replace('__OUTPUT_STEM__', outputStem)

    """
    set -euo pipefail
    ${prepareTss}

    computeMatrix reference-point \
        --referencePoint TSS \
        --beforeRegionStartLength 3000 \
        --afterRegionStartLength 3000 \
        --binSize 10 \
        --regionsFileName "${outputStem}.tss.bed" \
        --scoreFileName "coverage.RPKM.bw" \
        --numberOfProcessors "${task.cpus}" \
        --outFileName "${outputStem}.tss_matrix.gz" \
        --outFileNameMatrix "${outputStem}.tss_matrix.tsv"

    plotProfile \
        --matrixFile "${outputStem}.tss_matrix.gz" \
        --outFileName "${outputStem}.tss_profile.png" \
        --outFileNameData "${outputStem}.tss_profile.tsv" \
        --plotTitle "TSS enrichment"

    printf 'sample_id\\tannotation_mode\\tstatus\\n' \
        > "${outputStem}.tss_status.tsv"
    printf '%s\\t%s\\tcomputed\\n' "${sampleId}" "${annotation_mode}" \
        >> "${outputStem}.tss_status.tsv"

    printf 'TSS_ENRICHMENT:\\n  deeptools: ' \
        > "tss_enrichment_versions.yml"
    computeMatrix --version 2>&1 | sed -n '1p' \
        >> "tss_enrichment_versions.yml"
    """
}
