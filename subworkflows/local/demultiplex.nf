include { VALIDATE_MANIFEST } from '../../modules/local/validate_manifest'
include { DEMULTIPLEX_I2 } from '../../modules/local/demultiplex_i2'
include { FASTQC } from '../../modules/local/fastqc'

workflow DEMULTIPLEX {
    take:
    manifest
    barcode_mismatches
    allow_empty

    main:
    VALIDATE_MANIFEST(manifest)

    manifest_rows = VALIDATE_MANIFEST.out.normalized.flatMap { normalized_json ->
        def records = new groovy.json.JsonSlurper().parse(normalized_json.toFile())
        records.collect { raw_record ->
            def record = raw_record.collectEntries { key, value ->
                [(key.toString()): value]
            }
            tuple(record.library_id.toString(), record)
        }
    }

    library_inputs = manifest_rows
        .groupTuple(by: 0)
        .map { library_id, records ->
            def first = records.first()
            def library_meta = [library_id: library_id]
            tuple(
                library_meta,
                file(first.r1),
                file(first.r2),
                file(first.i2),
                records
            )
        }

    DEMULTIPLEX_I2(library_inputs, barcode_mismatches, allow_empty)

    derived_reads = DEMULTIPLEX_I2.out.reads.flatMap {
        library_meta, records, r1_files, r2_files ->
        def r1_list = r1_files instanceof java.util.List ? r1_files : [r1_files]
        def r2_list = r2_files instanceof java.util.List ? r2_files : [r2_files]

        records.collect { record ->
            def sample_id = record.sample_id.toString()
            def r1 = r1_list.find { path -> path.name == "${sample_id}_R1.fastq.gz" }
            def r2 = r2_list.find { path -> path.name == "${sample_id}_R2.fastq.gz" }
            if (r1 == null || r2 == null) {
                throw new IllegalStateException(
                    "demultiplexed FASTQs are missing for sample ${sample_id} " +
                    "from library ${library_meta.library_id}"
                )
            }
            def meta = new LinkedHashMap(record)
            meta.remove('r1')
            meta.remove('r2')
            meta.remove('i2')
            tuple(meta, r1, r2)
        }
    }

    FASTQC(derived_reads)

    versions_ch = VALIDATE_MANIFEST.out.versions.mix(
        DEMULTIPLEX_I2.out.versions,
        FASTQC.out.versions
    )

    emit:
    reads = derived_reads
    metrics = DEMULTIPLEX_I2.out.metrics
    fastqc = FASTQC.out.reports
    versions = versions_ch
}
