def validateNarrowMacsGenomeSize(rawValue) {
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

process MACS2_NARROW {
    tag "${meta.sample_id} vs ${control_meta.sample_id} (motif QC)"
    label 'process_heavy'

    conda "${projectDir}/envs/macs2.yml"
    container 'quay.io/biocontainers/macs2:2.2.9.1--py39hbcbf7aa_4'

    publishDir "${params.outdir}/peaks/${meta.sample_id}/narrow_motif_qc",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(target_bam, stageAs: 'target.analysis.bam'),
        path(target_bai, stageAs: 'target.analysis.bam.bai'),
        val(control_meta),
        path(control_bam, stageAs: 'control.analysis.bam'),
        path(control_bai, stageAs: 'control.analysis.bam.bai')
    val macs_genome_size
    val narrow_enabled

    output:
    tuple val(meta), path("sample_peaks.narrowPeak"), emit: motif_peaks
    tuple val(meta), path("sample_summits.bed"), emit: motif_summits
    tuple val(meta), path("sample_peaks.xls"), emit: auxiliary
    tuple val(meta), path("macs2_narrow.log"), emit: logs
    tuple val(meta), path("macs2_narrow_versions.yml"), emit: versions

    when:
    narrow_enabled == true

    script:
    def genomeSize = validateNarrowMacsGenomeSize(macs_genome_size)

    """
    set -euo pipefail
    macs2 callpeak \
        -t "target.analysis.bam" \
        -c "control.analysis.bam" \
        -g "${genomeSize}" \
        -f BAMPE \
        -n "sample" \
        --llocal 100000 \
        --keep-dup 1 \
        > "macs2_narrow.log" 2>&1

    touch "sample_peaks.narrowPeak"
    touch "sample_summits.bed"
    touch "sample_peaks.xls"

    printf 'MACS2_NARROW:\\n  macs2: ' > "macs2_narrow_versions.yml"
    macs2 --version 2>&1 | sed -n '1p' >> "macs2_narrow_versions.yml"
    """
}
