process BOWTIE2_BUILD {
    tag "${reference_meta.reference_id}"
    label 'process_heavy'

    conda "${projectDir}/envs/bowtie2.yml"
    container 'quay.io/biocontainers/bowtie2:2.5.4--he20e202_2'

    publishDir "${params.outdir}/pipeline_info/bowtie2_index",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename.tokenize('/').last() }

    input:
    tuple val(reference_meta), val(reference_mode), path(reference_inputs, stageAs: 'reference_input/*')

    output:
    tuple val(reference_meta), path("bt2_index/reference*.bt2*"), emit: index
    tuple val(reference_meta), path("bowtie2_build_versions.yml"), emit: versions

    script:
    if (!(reference_mode in ['fasta', 'index'])) {
        throw new IllegalArgumentException(
            "reference_mode must be 'fasta' or 'index', got ${reference_mode}"
        )
    }

    def prepareCommand
    if (reference_mode == 'fasta') {
        prepareCommand = '''
        reference_count=$(find "reference_input" -mindepth 1 -maxdepth 1 \
            ! -type d | wc -l | tr -d '[:space:]')
        if [[ "${reference_count}" -ne 1 ]]; then
            printf 'FASTA reference mode requires exactly one input; found %s\n' \
                "${reference_count}" >&2
            exit 1
        fi
        reference_file=$(find "reference_input" -mindepth 1 -maxdepth 1 \
            ! -type d -print -quit)
        bowtie2-build \
            --threads "${threads}" \
            "${reference_file}" \
            "bt2_index/reference"
        '''
    } else {
        prepareCommand = '''
        while IFS= read -r -d '' index_file; do
            index_name=$(basename "${index_file}")
            case "${index_name}" in
                *.rev.1.bt2)  suffix='rev.1.bt2' ;;
                *.rev.2.bt2)  suffix='rev.2.bt2' ;;
                *.1.bt2)      suffix='1.bt2' ;;
                *.2.bt2)      suffix='2.bt2' ;;
                *.3.bt2)      suffix='3.bt2' ;;
                *.4.bt2)      suffix='4.bt2' ;;
                *.rev.1.bt2l) suffix='rev.1.bt2l' ;;
                *.rev.2.bt2l) suffix='rev.2.bt2l' ;;
                *.1.bt2l)     suffix='1.bt2l' ;;
                *.2.bt2l)     suffix='2.bt2l' ;;
                *.3.bt2l)     suffix='3.bt2l' ;;
                *.4.bt2l)     suffix='4.bt2l' ;;
                *)
                    printf 'Unexpected Bowtie2 index component: %s\n' \
                        "${index_name}" >&2
                    exit 1
                    ;;
            esac
            destination="bt2_index/reference.${suffix}"
            if [[ -e "${destination}" ]]; then
                printf 'Duplicate Bowtie2 index component for %s\n' \
                    "${suffix}" >&2
                exit 1
            fi
            cp -- "${index_file}" "${destination}"
        done < <(find "reference_input" -mindepth 1 -maxdepth 1 \
            ! -type d -print0)

        small_count=$(find "bt2_index" -maxdepth 1 -type f -name 'reference*.bt2' | wc -l | tr -d '[:space:]')
        large_count=$(find "bt2_index" -maxdepth 1 -type f -name 'reference*.bt2l' | wc -l | tr -d '[:space:]')
        if [[ "${small_count}" -eq 6 && "${large_count}" -eq 0 ]]; then
            index_extension='bt2'
        elif [[ "${large_count}" -eq 6 && "${small_count}" -eq 0 ]]; then
            index_extension='bt2l'
        else
            printf 'Supplied Bowtie2 index must contain exactly one complete six-file .bt2 or .bt2l set; found %s small and %s large components\n' \
                "${small_count}" "${large_count}" >&2
            exit 1
        fi

        for index_suffix in 1 2 3 4 rev.1 rev.2; do
            if [[ ! -f "bt2_index/reference.${index_suffix}.${index_extension}" ]]; then
                printf 'Missing Bowtie2 index component: %s\n' \
                    "bt2_index/reference.${index_suffix}.${index_extension}" >&2
                exit 1
            fi
        done
        '''
    }

    """
    set -euo pipefail
    threads="${task.cpus}"
    mkdir -p "bt2_index"
    ${prepareCommand}

    printf 'BOWTIE2_BUILD:\\n  bowtie2: ' > "bowtie2_build_versions.yml"
    bowtie2 --version 2>&1 | sed -n '1p' >> "bowtie2_build_versions.yml"
    """
}
