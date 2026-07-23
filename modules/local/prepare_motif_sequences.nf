def motifSequenceSampleId(meta) {
    def safeId = /[A-Za-z0-9][A-Za-z0-9._-]*/
    if (!(meta instanceof Map)) {
        throw new IllegalArgumentException(
            "motif-sequence metadata must be a map"
        )
    }
    def sampleId = meta.sample_id?.toString()
    if (sampleId == null || !sampleId.matches(safeId)) {
        throw new IllegalArgumentException(
            "motif-sequence sample_id must match ${safeId}, got ${sampleId}"
        )
    }
    sampleId
}

def validateMotifWindow(rawWindow) {
    def text = rawWindow?.toString()
    if (text == null || !(text ==~ /[0-9]+/)) {
        throw new IllegalArgumentException(
            "motif_window must be a positive integer, got ${text}"
        )
    }
    def value = new BigInteger(text)
    if (value.signum() <= 0 || value > Integer.MAX_VALUE) {
        throw new IllegalArgumentException(
            "motif_window must be a positive integer, got ${text}"
        )
    }
    value.toString()
}

def validateMotifAnchorMode(rawMode) {
    def mode = rawMode?.toString()
    if (!(mode in ['summit', 'midpoint'])) {
        throw new IllegalArgumentException(
            "unsupported motif anchor mode ${mode}; expected summit or midpoint"
        )
    }
    mode
}

process PREPARE_MOTIF_SEQUENCES {
    tag "${meta.sample_id} (${anchor_mode})"
    label 'process_standard'

    conda "${projectDir}/envs/bedtools.yml"
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'

    publishDir "${params.outdir}/motifs/${meta.sample_id}/sequences",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta),
        path(anchor_regions, stageAs: 'motif.anchor.bed'),
        val(anchor_mode),
        path(final_peaks, stageAs: 'final.peaks.bed')
    path fasta, stageAs: 'reference.fa'
    path(blacklist_files, stageAs: 'blacklist/regions*.bed', arity: '0..1')
    val motif_window

    output:
    tuple val(meta),
        path("foreground.fa"),
        path("background.fa"),
        path("motif_sequence_status.tsv"),
        emit: sequences
    tuple val(meta),
        path("foreground.windows.bed"),
        path("background.windows.bed"),
        emit: windows
    tuple val(meta), path("prepare_motif_sequences.log"), emit: logs
    tuple val(meta), path("prepare_motif_sequences_versions.yml"), emit: versions

    script:
    motifSequenceSampleId(meta)
    def window = validateMotifWindow(motif_window)
    def anchorMode = validateMotifAnchorMode(anchor_mode)
    def suppliedBlacklists
    if (blacklist_files == null) {
        suppliedBlacklists = []
    } else if (blacklist_files instanceof java.util.Collection) {
        suppliedBlacklists = blacklist_files as List
    } else {
        suppliedBlacklists = [blacklist_files]
    }
    if (suppliedBlacklists.size() > 1) {
        throw new IllegalArgumentException(
            "PREPARE_MOTIF_SEQUENCES accepts at most one blacklist, got " +
            suppliedBlacklists.size()
        )
    }
    def blacklistAppend = suppliedBlacklists
        ? 'cat "blacklist/regions.bed" >> "background_exclusions.unsorted.bed"'
        : 'true'
    def anchorPreparation
    if (anchorMode == 'summit') {
        def blacklistSummitFilter = suppliedBlacklists
            ? '''
    bedtools intersect -v -a "anchors.in.final.peaks.bed" -b "blacklist/regions.bed" \
        > "usable.anchors.bed" 2>> "prepare_motif_sequences.log"
    '''
            : 'cp "anchors.in.final.peaks.bed" "usable.anchors.bed"'
        anchorPreparation = """
    bedtools intersect -u -a "motif.anchor.bed" -b "final.peaks.bed" \
        > "anchors.in.final.peaks.bed" 2>> "prepare_motif_sequences.log"
    ${blacklistSummitFilter}
    """
    } else {
        anchorPreparation = 'cp "final.peaks.bed" "usable.anchors.bed"'
    }

    """
    set -euo pipefail
    : > "prepare_motif_sequences.log"

    awk '
        BEGIN { OFS = "\\t"; name = ""; sequence_length = 0 }
        /^>/ {
            if (name != "") {
                print name, sequence_length
            }
            header = substr(\$0, 2)
            sub(/[[:space:]].*$/, "", header)
            if (header == "") {
                print "reference FASTA contains a blank sequence name" > "/dev/stderr"
                exit 2
            }
            name = header
            sequence_length = 0
            next
        }
        {
            gsub(/[[:space:]]/, "", \$0)
            sequence_length += length(\$0)
        }
        END {
            if (name != "") {
                print name, sequence_length
            }
        }
    ' "reference.fa" > "reference.genome"

    if [[ ! -s "reference.genome" ]]; then
        printf '%s\\n' 'reference FASTA contains no sequences' >&2
        exit 2
    fi

    ${anchorPreparation}

    cat > "make_windows.awk" <<'AWK'
BEGIN {
    OFS = "\t"
    count = 0
}
FNR == NR {
    if (NF != 2 || $2 !~ /^[0-9]+$/ || $2 <= 0) {
        print "invalid reference genome-size row " FNR > "/dev/stderr"
        exit 2
    }
    chromosome_size[$1] = $2
    next
}
/^[[:space:]]*($|#)/ {
    next
}
{
    if (NF < 3 || !($1 in chromosome_size)) {
        print "motif anchor row " FNR " has an unknown chromosome" > "/dev/stderr"
        exit 2
    }
    if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ ||
        $2 < 0 || $3 <= $2 || $3 > chromosome_size[$1]) {
        print "motif anchor row " FNR " has invalid coordinates" > "/dev/stderr"
        exit 2
    }
    center = int(($2 + $3) / 2)
    requested_start = center - int(window / 2)
    maximum_start = chromosome_size[$1] - window
    if (maximum_start < 0) {
        start = 0
        end = chromosome_size[$1]
    } else {
        start = requested_start
        if (start < 0) {
            start = 0
        }
        if (start > maximum_start) {
            start = maximum_start
        }
        end = start + window
    }
    count += 1
    printf "%s\t%d\t%d\tpeak_%06d\n", $1, start, end, count
}
AWK

    awk \
        -v mode="${anchorMode}" \
        -v window="${window}" \
        -f "make_windows.awk" \
        "reference.genome" \
        "usable.anchors.bed" \
        > "foreground.windows.bed"

    : > "foreground.fa"
    : > "background.windows.bed"
    : > "background.fa"

    if [[ ! -s "foreground.windows.bed" ]]; then
        printf 'metric\\tvalue\\nstatus\\tno_peaks\\nanchor_mode\\t%s\\nwindow\\t%s\\n' \
            "${anchorMode}" "${window}" \
            > "motif_sequence_status.tsv"
        printf '%s\\n' \
            'no usable peak anchors; motif sequence preparation skipped' \
            >> "prepare_motif_sequences.log"
    else
        bedtools getfasta \
            -fi "reference.fa" \
            -bed "foreground.windows.bed" \
            -nameOnly \
            > "foreground.fa" \
            2>> "prepare_motif_sequences.log"

        cp "final.peaks.bed" "background_exclusions.unsorted.bed"
        ${blacklistAppend}
        bedtools sort \
            -i "background_exclusions.unsorted.bed" \
            2>> "prepare_motif_sequences.log" |
        bedtools merge \
            -i - \
            > "background_exclusions.bed" \
            2>> "prepare_motif_sequences.log"

        awk '
            BEGIN { OFS = "\\t"; copies = 8; candidate = 0 }
            {
                for (copy = 1; copy <= copies; copy++) {
                    candidate += 1
                    print \$1, \$2, \$3, "candidate_" candidate
                }
            }
        ' "foreground.windows.bed" > "candidate.templates.bed"

        bedtools shuffle \
            -i "candidate.templates.bed" \
            -g "reference.genome" \
            -excl "background_exclusions.bed" \
            -seed 1729 \
            -noOverlapping \
            -maxTries 10000 \
            > "candidate.windows.bed" \
            2>> "prepare_motif_sequences.log"

        bedtools getfasta \
            -fi "reference.fa" \
            -bed "foreground.windows.bed" \
            -nameOnly \
            -tab \
            > "foreground.sequences.tsv" \
            2>> "prepare_motif_sequences.log"
        bedtools getfasta \
            -fi "reference.fa" \
            -bed "candidate.windows.bed" \
            -nameOnly \
            -tab \
            > "candidate.sequences.tsv" \
            2>> "prepare_motif_sequences.log"

        awk '
            BEGIN { OFS = "\\t" }
            {
                sequence = toupper(\$2)
                sequence_length = length(sequence)
                gc_count = gsub(/[GC]/, "", sequence)
                if (sequence_length == 0) {
                    print "empty foreground sequence for " \$1 > "/dev/stderr"
                    exit 2
                }
                print \$1, sequence_length, gc_count / sequence_length
            }
        ' "foreground.sequences.tsv" > "foreground.features.tsv"
        awk '
            BEGIN { OFS = "\\t" }
            {
                sequence = toupper(\$2)
                sequence_length = length(sequence)
                gc_count = gsub(/[GC]/, "", sequence)
                if (sequence_length == 0) {
                    print "empty background candidate for " \$1 > "/dev/stderr"
                    exit 2
                }
                print \$1, sequence_length, gc_count / sequence_length
            }
        ' "candidate.sequences.tsv" > "candidate.sequence_features.tsv"
        paste "candidate.windows.bed" "candidate.sequence_features.tsv" |
        awk '
            BEGIN { OFS = "\\t" }
            {
                if (\$4 != \$5) {
                    print "candidate BED/FASTA order mismatch" > "/dev/stderr"
                    exit 2
                }
                print \$1, \$2, \$3, \$4, \$6, \$7
            }
        ' > "candidate.features.tsv"

        cat > "match_gc.awk" <<'AWK'
BEGIN {
    OFS = "\t"
    foreground_count = 0
    candidate_count = 0
}
FNR == NR {
    foreground_count += 1
    foreground_name[foreground_count] = $1
    foreground_length[foreground_count] = $2
    foreground_gc[foreground_count] = $3
    next
}
{
    candidate_count += 1
    candidate_chromosome[candidate_count] = $1
    candidate_start[candidate_count] = $2
    candidate_end[candidate_count] = $3
    candidate_name[candidate_count] = $4
    candidate_length[candidate_count] = $5
    candidate_gc[candidate_count] = $6
}
END {
    for (foreground_index = 1;
         foreground_index <= foreground_count;
         foreground_index++) {
        best = 0
        best_distance = 2
        for (candidate_index = 1;
             candidate_index <= candidate_count;
             candidate_index++) {
            if (used[candidate_index] ||
                candidate_length[candidate_index] != foreground_length[foreground_index]) {
                continue
            }
            distance = candidate_gc[candidate_index] - foreground_gc[foreground_index]
            if (distance < 0) {
                distance = -distance
            }
            if (best == 0 || distance < best_distance) {
                best = candidate_index
                best_distance = distance
            }
        }
        if (best == 0) {
            print "no length-matched shuffled background for " \
                foreground_name[foreground_index] > "/dev/stderr"
            exit 2
        }
        used[best] = 1
        print candidate_chromosome[best], candidate_start[best], \
            candidate_end[best], foreground_name[foreground_index]
    }
}
AWK

        awk \
            -f "match_gc.awk" \
            "foreground.features.tsv" \
            "candidate.features.tsv" \
            > "background.windows.bed"

        bedtools getfasta \
            -fi "reference.fa" \
            -bed "background.windows.bed" \
            -nameOnly \
            > "background.fa" \
            2>> "prepare_motif_sequences.log"

        printf 'metric\\tvalue\\nstatus\\tready\\nanchor_mode\\t%s\\nwindow\\t%s\\nseed\\t1729\\n' \
            "${anchorMode}" "${window}" \
            > "motif_sequence_status.tsv"
    fi

    printf 'PREPARE_MOTIF_SEQUENCES:\\n  bedtools: ' \
        > "prepare_motif_sequences_versions.yml"
    bedtools --version 2>&1 | sed -n '1p' \
        >> "prepare_motif_sequences_versions.yml"
    """
}
