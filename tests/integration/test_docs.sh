#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${repo_root}"

python3 - <<'PY'
import json
import re
import sys
from pathlib import Path

root = Path(".")
readme_path = root / "README.md"
if not readme_path.is_file():
    raise SystemExit("FAIL: missing README.md")

readme = readme_path.read_text(encoding="utf-8")
config = (root / "nextflow.config").read_text(encoding="utf-8")
schema = json.loads((root / "nextflow_schema.json").read_text(encoding="utf-8"))
main = (root / "main.nf").read_text(encoding="utf-8")
qc = (root / "subworkflows/local/qc.nf").read_text(encoding="utf-8")
multiqc = (root / "modules/local/multiqc.nf").read_text(encoding="utf-8")

minimum_nextflow = "23.10.0"
version_declarations = re.findall(
    r"(?m)^\s*manifest\.nextflowVersion\s*=\s*'([^']+)'\s*$",
    config,
)
if version_declarations != [f">={minimum_nextflow}"]:
    raise SystemExit(
        "FAIL: nextflow.config must enforce exactly "
        f"manifest.nextflowVersion = '>={minimum_nextflow}'"
    )
minimum_version_phrase = f"Nextflow {minimum_nextflow} or later"
if minimum_version_phrase not in readme:
    raise SystemExit(
        "FAIL: README must document the exact minimum runtime as "
        + minimum_version_phrase
    )

workflow_text = "\n".join(
    path.read_text(encoding="utf-8")
    for source_root in (root / "modules", root / "subworkflows")
    for path in sorted(source_root.rglob("*.nf"))
)
invoked_bin_tools = sorted(
    set(
        re.findall(
            r"(?:\$\{projectDir\}/bin/|\$\{projectBin\}/|bin/)"
            r"([A-Za-z0-9_.-]+\.py)\b",
            workflow_text,
        )
    )
)
if not invoked_bin_tools:
    raise SystemExit("FAIL: no repository Python tools were discovered in workflow sources")
missing_bin_tools = [
    tool for tool in invoked_bin_tools if not (root / "bin" / tool).is_file()
]
if missing_bin_tools:
    raise SystemExit(
        "FAIL: workflow sources invoke missing repository bin tools: "
        + ", ".join(missing_bin_tools)
    )

params_block = re.search(r"(?ms)^\s*params\s*\{(.*?)^\s*\}", config)
if params_block is None:
    raise SystemExit("FAIL: nextflow.config has no params block")
config_params = set(
    re.findall(r"(?m)^\s{4}([A-Za-z_][A-Za-z0-9_]*)\s*=", params_block.group(1))
)
schema_params = set(schema.get("properties", {}))
if config_params != schema_params:
    missing_schema = sorted(config_params - schema_params)
    missing_config = sorted(schema_params - config_params)
    raise SystemExit(
        "FAIL: config/schema parameter drift; "
        f"missing from schema={missing_schema}, missing from config={missing_config}"
    )

missing_params = sorted(
    name for name in config_params if f"--{name}" not in readme
)
if missing_params:
    raise SystemExit(
        "FAIL: README does not document CLI parameters: "
        + ", ".join(f"--{name}" for name in missing_params)
    )

environment_files = sorted(path.relative_to(root).as_posix() for path in (root / "envs").glob("*.yml"))
missing_envs = [path for path in environment_files if path not in readme]
if missing_envs:
    raise SystemExit(
        "FAIL: README does not document environment files: " + ", ".join(missing_envs)
    )

required_manifest_columns = [
    "sample_id",
    "library_id",
    "input_group",
    "barcode",
    "assay_target",
    "is_control",
    "control_id",
    "expected_motif",
    "r1",
    "r2",
    "i2",
]
missing_columns = [
    column for column in required_manifest_columns
    if f"`{column}`" not in readme
]
if missing_columns:
    raise SystemExit(
        "FAIL: README does not document manifest columns: " + ", ".join(missing_columns)
    )

required_environment_variables = ["NXF_CONDA_CACHEDIR", "NXF_HOME"]
missing_variables = [
    name for name in required_environment_variables if name not in readme
]
if missing_variables:
    raise SystemExit(
        "FAIL: README does not document runtime environment variables: "
        + ", ".join(missing_variables)
    )

required_phrases = [
    "TATAGCCT",
    "ATAGAGGC",
    "NX701",
    "NX702",
    "NX703",
    "NX704",
    "NX705",
    "NX706",
    "FRiP",
    "duplicate",
    "MAPQ",
    "TSS enrichment",
    "broad",
    "narrow",
    "JASPAR",
    "-resume",
    "Nextflow runtime execution was not verified in this workspace",
]
missing_phrases = [phrase for phrase in required_phrases if phrase not in readme]
if missing_phrases:
    raise SystemExit(
        "FAIL: README is missing required documentation topics: "
        + ", ".join(missing_phrases)
    )

dashboard_checks = {
    "QC collects optional enrichment files and forwards them to MultiQC":
        "enrichment_files" in qc
        and "safe_enrichment_files = enrichment_files" in qc
        and ".ifEmpty { ignored -> [] }" in qc
        and re.search(
            r"MULTIQC\(\s*fastqc_files,\s*demux_custom_files,\s*"
            r"library_custom_files,\s*insert_size_files,\s*peak_qc_files,\s*"
            r"motif_metric_files,\s*tss_status_files,\s*safe_enrichment_files,\s*"
            r"annotation_status\s*\)",
            qc,
            re.S,
        ),
    "MultiQC stages optional enrichment inputs, publishes enrichment artifacts, and skips empty custom sections":
        "path enrichment_files, stageAs: 'enrichment??/*'" in multiqc
        and "nanocut_peak_enrichment_mqc.tsv" in multiqc
        and "if enrichment_rows:" in multiqc
        and "copyfile(" in multiqc
        and "publishDir" in multiqc
        and "pattern: 'peak_enrichment.tsv'" in multiqc
        and "pattern: 'observed_vs_null.png'" in multiqc,
    "MultiQC adds a dedicated peak-enrichment section with stable heatmap references for all four null models":
        'section_name: "Nano-CUT&Tag peak enrichment"' in multiqc
        and all(model in multiqc for model in (
            "random",
            "length_matched",
            "gc_matched",
            "length_gc_matched",
        ))
        and all(name in multiqc for name in (
            "matrix_random.png",
            "matrix_length_matched.png",
            "matrix_gc_matched.png",
            "matrix_length_gc_matched.png",
        )),
    "Completion summary links enrichment TSV only when the dashboard artifacts exist":
        "path enrichment_files, stageAs: 'enrichment??/*'" in main
        and "reports/summary/peak_enrichment.tsv" in main
        and "if enrichment_tables:" in main,
}

dashboard_failures = [
    description for description, passed in dashboard_checks.items() if not passed
]
if dashboard_failures:
    raise SystemExit(
        "FAIL: dashboard structural assertions failed:\n  - "
        + "\n  - ".join(dashboard_failures)
    )

print(
    f"README documents all {len(config_params)} CLI parameters, "
    f"{len(environment_files)} Conda environments, manifest fields, "
    "runtime environment variables, and required operational/QC topics; "
    f"all {len(invoked_bin_tools)} referenced repository bin tools exist."
)
PY
