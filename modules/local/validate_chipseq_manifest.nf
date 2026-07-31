process VALIDATE_CHIPSEQ_MANIFEST {
    tag "${manifest.simpleName}"
    label 'process_light'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy',
        overwrite: true

    input:
    path manifest, stageAs: 'chipseq_manifest.csv'
    path fasta, stageAs: 'reference.fa'

    output:
    path "normalized_chipseq_manifest.tsv", emit: normalized
    path "validate_chipseq_manifest_versions.yml", emit: versions

    script:
    def projectBin = "${projectDir}/bin".toString()

    """
    set -euo pipefail
    python <<'PY'
import csv
import sys
from pathlib import Path

sys.path.insert(0, r"${projectBin}")

from peak_enrichment import load_fasta_sequences, load_reference_manifest

fasta_sequences = load_fasta_sequences("reference.fa")
chrom_sizes = {
    chrom: len(sequence)
    for chrom, sequence in fasta_sequences.items()
}
records = load_reference_manifest("chipseq_manifest.csv", chrom_sizes)

with open("normalized_chipseq_manifest.tsv", "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["reference_id", "tf", "reference_type", "peak_file"],
        delimiter="\\t",
    )
    writer.writeheader()
    for record in records:
        writer.writerow(
            {
                "reference_id": record.reference_id,
                "tf": record.tf,
                "reference_type": record.reference_type,
                "peak_file": str(record.peak_file),
            }
        )
PY

    printf 'VALIDATE_CHIPSEQ_MANIFEST:\\n  python: ' \
        > "validate_chipseq_manifest_versions.yml"
    python --version 2>&1 >> "validate_chipseq_manifest_versions.yml"
    printf '  peak_enrichment.py: repository\\n' \
        >> "validate_chipseq_manifest_versions.yml"
    """
}
