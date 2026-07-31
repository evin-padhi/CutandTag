process VALIDATE_CHIPSEQ_MANIFEST {
    tag 'chipseq-manifest'
    label 'process_light'

    conda "${projectDir}/envs/python.yml"
    container 'python:3.12.3-slim-bookworm'

    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy',
        overwrite: true,
        pattern: 'normalized_chipseq_manifest.tsv'
    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy',
        overwrite: true,
        pattern: 'validate_chipseq_manifest_versions.yml'

    input:
    tuple val(reference_rows_b64),
        path(reference_peak_files, stageAs: 'incoming_reference??????/*', arity: '1..*')
    path fasta, stageAs: 'reference.fa'
    path enrichment_script, stageAs: 'bin/peak_enrichment.py'

    output:
    path "normalized_chipseq_manifest.tsv", emit: normalized
    path "reference_peaks/*.bed", emit: peaks
    path "validate_chipseq_manifest_versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    printf '%s' "${reference_rows_b64}" | python -c \
        'import base64, sys; sys.stdout.buffer.write(base64.b64decode(sys.stdin.buffer.read()))' \
        > "reference_rows.json"

    python <<'PY'
import csv
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path("bin").resolve()))

from peak_enrichment import load_fasta_chrom_sizes, load_public_chipseq_manifest

chrom_sizes = load_fasta_chrom_sizes("reference.fa")
rows = json.loads(Path("reference_rows.json").read_text(encoding="utf-8"))
staged_peaks = sorted(Path().glob("incoming_reference*/*"))
if len(rows) != len(staged_peaks):
    raise SystemExit(
        "reference metadata/file count mismatch: "
        f"{len(rows)} rows for {len(staged_peaks)} BED files"
    )

with open("staged_chipseq_manifest.csv", "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["reference_id", "tf", "peak_file"],
    )
    writer.writeheader()
    for row, peak_file in zip(rows, staged_peaks, strict=True):
        writer.writerow(
            {
                "reference_id": row["reference_id"],
                "tf": row["tf"],
                "peak_file": str(peak_file),
            }
        )

records = load_public_chipseq_manifest("staged_chipseq_manifest.csv", chrom_sizes)
output_dir = Path("reference_peaks")
output_dir.mkdir()

with open("normalized_chipseq_manifest.tsv", "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["reference_id", "tf", "reference_type", "peak_file"],
        delimiter="\\t",
    )
    writer.writeheader()
    for record in records:
        destination = output_dir / f"{record.reference_id}.bed"
        shutil.copyfile(record.peak_file, destination)
        writer.writerow(
            {
                "reference_id": record.reference_id,
                "tf": record.tf,
                "reference_type": "chipseq",
                "peak_file": f"reference_peaks/{destination.name}",
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
