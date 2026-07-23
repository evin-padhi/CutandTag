process VALIDATE_MANIFEST {
    tag "${manifest.simpleName}"
    label 'process_light'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy',
        overwrite: true

    input:
    path manifest, stageAs: 'input_manifest.csv'

    output:
    path "normalized_manifest.json", emit: normalized
    path "validate_manifest_versions.yml", emit: versions

    script:
    """
    manifest.py validate \
        --input "input_manifest.csv" \
        --output normalized_manifest.json

    printf 'VALIDATE_MANIFEST:\\n  python: ' > validate_manifest_versions.yml
    python --version 2>&1 >> validate_manifest_versions.yml
    printf '  manifest.py: repository\\n' >> validate_manifest_versions.yml
    """
}
