include { BOWTIE2_BUILD } from '../../modules/local/bowtie2_build'
include { BOWTIE2_ALIGN } from '../../modules/local/bowtie2_align'
include { SAMTOOLS_SORT_INDEX } from '../../modules/local/samtools_sort_index'
include { SAMTOOLS_FILTER } from '../../modules/local/samtools_filter'
include { SAMTOOLS_METRICS } from '../../modules/local/samtools_metrics'
include { BAMCOVERAGE } from '../../modules/local/bamcoverage'

def validateAlignmentMeta(meta) {
    def SAFE_ID = /[A-Za-z0-9][A-Za-z0-9._-]*/
    if (!(meta instanceof Map)) {
        throw new IllegalArgumentException(
            "alignment metadata must be a map, got ${meta?.getClass()?.name}"
        )
    }
    if (!meta.containsKey('sample_id')) {
        throw new IllegalArgumentException(
            "alignment metadata is missing required key sample_id"
        )
    }
    def sampleId = meta.sample_id?.toString()
    if (sampleId == null || !sampleId.matches(SAFE_ID)) {
        throw new IllegalArgumentException(
            "sample_id must match ${SAFE_ID}, got ${sampleId}"
        )
    }
    new LinkedHashMap(meta)
}

workflow ALIGN_QC {
    take:
    reads
    fasta
    bowtie2_index
    min_mapq

    main:
    safe_reads = reads.map { meta, r1, r2 ->
        tuple(validateAlignmentMeta(meta), r1, r2)
    }

    reference_meta = [reference_id: 'bowtie2']
    if (bowtie2_index) {
        index_pattern = "${bowtie2_index}.{1,2,3,4,rev.1,rev.2}.bt2{,l}"
        reference_source = Channel
            .fromPath(index_pattern, checkIfExists: true)
            .collect()
            .map { index_files ->
                tuple(reference_meta, 'index', index_files)
            }
        reference_mode = 'index'
    } else if (fasta) {
        reference_source = Channel
            .fromPath(fasta, checkIfExists: true)
            .map { reference_fasta ->
                tuple(reference_meta, 'fasta', [reference_fasta])
            }
        reference_mode = 'fasta'
    } else {
        throw new IllegalArgumentException(
            "ALIGN_QC requires either fasta or bowtie2_index"
        )
    }

    BOWTIE2_BUILD(reference_source)

    alignment_inputs = safe_reads.combine(BOWTIE2_BUILD.out.index)
    BOWTIE2_ALIGN(alignment_inputs)
    SAMTOOLS_SORT_INDEX(BOWTIE2_ALIGN.out.sam)
    SAMTOOLS_FILTER(SAMTOOLS_SORT_INDEX.out.bam, min_mapq)
    SAMTOOLS_METRICS(SAMTOOLS_SORT_INDEX.out.bam)
    BAMCOVERAGE(SAMTOOLS_FILTER.out.bam, min_mapq)

    versions_ch = BOWTIE2_BUILD.out.versions.mix(
        BOWTIE2_ALIGN.out.versions,
        SAMTOOLS_SORT_INDEX.out.versions,
        SAMTOOLS_FILTER.out.versions,
        SAMTOOLS_METRICS.out.versions,
        BAMCOVERAGE.out.versions
    )

    emit:
    index = BOWTIE2_BUILD.out.index
    analysis_bam = SAMTOOLS_SORT_INDEX.out.bam
    filtered_bam = SAMTOOLS_FILTER.out.bam
    metrics = SAMTOOLS_METRICS.out.metrics
    coverage = BAMCOVERAGE.out.bigwig
    alignment_summary = BOWTIE2_ALIGN.out.summary
    versions = versions_ch
}
