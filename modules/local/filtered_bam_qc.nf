process FILTERED_BAM_QC {
    tag "${meta.sample_id}"
    label 'process_light'

    conda "${projectDir}/envs/samtools.yml"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    publishDir "${params.outdir}/qc/library/${meta.sample_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(filtered_bam, stageAs: 'qc.filtered.bam'),
        path(filtered_bai, stageAs: 'qc.filtered.bam.bai')

    output:
    tuple val(meta), path("filtered_bam_qc.tsv"), emit: metrics
    tuple val(meta), path("filtered_bam_qc_versions.yml"), emit: versions

    script:
    """
    set -euo pipefail
    cat > "filtered_flags.awk" <<'AWK'
function has_flag(flag, bit) {
    return int(flag / bit) % 2
}

{
    if (NF < 2 || \$2 !~ /^[0-9]+\$/) {
        print "filtered BAM contains a malformed SAM FLAG" > "/dev/stderr"
        exit 2
    }
    first_bit = has_flag(\$2, 64)
    second_bit = has_flag(\$2, 128)
    if (first_bit == second_bit) {
        print "each filtered alignment must set exactly one of 0x40 and 0x80" \
            " for " \$1 > "/dev/stderr"
        exit 2
    }
}
AWK

    samtools view "qc.filtered.bam" \
        | LC_ALL=C awk -f "filtered_flags.awk"

    filtered_reads=\$(samtools view -c "qc.filtered.bam")
    first_mates=\$(samtools view -c -f 64 -F 128 "qc.filtered.bam")
    second_mates=\$(samtools view -c -f 128 -F 64 "qc.filtered.bam")

    if [[ "\${first_mates}" -ne "\${second_mates}" ]]; then
        printf 'filtered first/second mate counts differ: first=%s second=%s\\n' \
            "\${first_mates}" "\${second_mates}" >&2
        exit 2
    fi
    if [[ \$((first_mates + second_mates)) -ne "\${filtered_reads}" ]]; then
        printf 'first-only plus second-only counts do not equal total records: first=%s second=%s total=%s\\n' \
            "\${first_mates}" "\${second_mates}" "\${filtered_reads}" >&2
        exit 2
    fi
    filtered_fragments="\${first_mates}"
    if [[ "\${filtered_reads}" -ne \$((filtered_fragments * 2)) ]]; then
        printf 'filtered read count is not twice the fragment count: reads=%s fragments=%s\\n' \
            "\${filtered_reads}" "\${filtered_fragments}" >&2
        exit 2
    fi

    printf 'metric\\tvalue\\n' > "filtered_bam_qc.tsv"
    printf 'filtered_reads\\t%s\\n' "\${filtered_reads}" \
        >> "filtered_bam_qc.tsv"
    printf 'filtered_fragments\\t%s\\n' "\${filtered_fragments}" \
        >> "filtered_bam_qc.tsv"

    printf 'FILTERED_BAM_QC:\\n  samtools: ' \
        > "filtered_bam_qc_versions.yml"
    samtools --version 2>&1 | sed -n '1p' \
        >> "filtered_bam_qc_versions.yml"
    """
}
