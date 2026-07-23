process BAM_TO_FRAGMENTS {
    tag "${meta.sample_id}"
    label 'process_standard'

    conda "${projectDir}/envs/samtools.yml"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    publishDir "${params.outdir}/qc/fragments/${meta.sample_id}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(filtered_bam, stageAs: 'qc.filtered.bam'),
        path(filtered_bai, stageAs: 'qc.filtered.bam.bai')

    output:
    tuple val(meta), path("fragments.bedpe"), emit: fragments
    tuple val(meta), path("bam_to_fragments_versions.yml"), emit: versions

    script:
    """
    set -euo pipefail

    cat > "fragment_pairs.awk" <<'AWK'
BEGIN {
    OFS = "\t"
    current = ""
    count = 0
    first_mate_count = 0
    second_mate_count = 0
    failed = 0
}

function die(message) {
    print "bam_to_fragments: " message > "/dev/stderr"
    failed = 1
    exit 2
}

function reference_end(position, cigar,    rest, token, length_value, operation, consumed) {
    if (cigar == "*" || position < 1) {
        die("mapped primary alignment has no usable CIGAR/POS for " current)
    }
    rest = cigar
    consumed = 0
    while (match(rest, /^[0-9]+[MIDNSHP=X]/)) {
        token = substr(rest, RSTART, RLENGTH)
        length_value = token + 0
        operation = substr(token, length(token), 1)
        if (index("MDN=X", operation) > 0) {
            consumed += length_value
        }
        rest = substr(rest, RLENGTH + 1)
    }
    if (rest != "" || consumed < 1) {
        die("invalid or zero-reference-width CIGAR for " current ": " cigar)
    }
    return position - 1 + consumed
}

function has_flag(flag, bit) {
    return int(flag / bit) % 2
}

function flush_pair(    swap, c1, s1, e1, c2, s2, e2) {
    if (count != 2) {
        die("expected exactly two primary proper-pair alignments for " current \
            ", observed " count)
    }
    if (first_mate_count != 1 || second_mate_count != 1) {
        die("expected one 0x40 and one 0x80 mate for " current \
            ", observed 0x40=" first_mate_count \
            " and 0x80=" second_mate_count)
    }
    swap = chromosome[2] < chromosome[1] \
        || (chromosome[2] == chromosome[1] && start[2] < start[1]) \
        || (chromosome[2] == chromosome[1] && start[2] == start[1] \
            && end[2] < end[1])
    if (swap) {
        c1 = chromosome[2]; s1 = start[2]; e1 = end[2]
        c2 = chromosome[1]; s2 = start[1]; e2 = end[1]
    } else {
        c1 = chromosome[1]; s1 = start[1]; e1 = end[1]
        c2 = chromosome[2]; s2 = start[2]; e2 = end[2]
    }
    print c1, s1, e1, c2, s2, e2, current
    count = 0
    first_mate_count = 0
    second_mate_count = 0
}

{
    if (NF < 6) {
        die("malformed SAM alignment")
    }
    if (current != "" && \$1 != current) {
        flush_pair()
    }
    if (current == "" || \$1 != current) {
        current = \$1
    }
    if (count >= 2) {
        die("more than two primary proper-pair alignments for " current)
    }
    if (\$2 !~ /^[0-9]+\$/) {
        die("SAM FLAG must be a non-negative integer for " current)
    }
    first_bit = has_flag(\$2, 64)
    second_bit = has_flag(\$2, 128)
    if (first_bit == second_bit) {
        die("each alignment must set exactly one of 0x40 and 0x80 for " current)
    }
    first_mate_count += first_bit
    second_mate_count += second_bit
    count += 1
    chromosome[count] = \$3
    if (chromosome[count] == "*") {
        die("proper-pair alignment is unexpectedly unmapped for " current)
    }
    start[count] = \$4 - 1
    end[count] = reference_end(\$4, \$6)
}

END {
    if (!failed && current != "") {
        flush_pair()
    }
}
AWK

    samtools sort \
        -n \
        -O BAM \
        -@ "${task.cpus}" \
        "qc.filtered.bam" \
        | samtools view \
            -f 2 \
            -F 2304 \
            - \
        > "name_sorted.primary_proper.sam"

    LC_ALL=C awk \
        -f "fragment_pairs.awk" \
        "name_sorted.primary_proper.sam" \
        | LC_ALL=C sort \
            -t \$'\\t' \
            -k1,1 \
            -k2,2n \
            -k3,3n \
            -k4,4 \
            -k5,5n \
            -k6,6n \
            -k7,7 \
        > "fragments.bedpe"

    printf 'BAM_TO_FRAGMENTS:\\n  samtools: ' \
        > "bam_to_fragments_versions.yml"
    samtools --version 2>&1 | sed -n '1p' \
        >> "bam_to_fragments_versions.yml"
    """
}
