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

module_text = "\n".join(
    path.read_text(encoding="utf-8")
    for path in sorted((root / "modules").rglob("*.nf"))
)
invoked_bin_tools = sorted(set(re.findall(r"\b([A-Za-z0-9_.-]+\.py)\b", module_text)))
if not invoked_bin_tools:
    raise SystemExit("FAIL: no repository Python tools were discovered in modules")
missing_bin_tools = [
    tool for tool in invoked_bin_tools if not (root / "bin" / tool).is_file()
]
if missing_bin_tools:
    raise SystemExit(
        "FAIL: modules invoke missing repository bin tools: "
        + ", ".join(missing_bin_tools)
    )
non_executable_bin_tools = [
    tool
    for tool in invoked_bin_tools
    if not ((root / "bin" / tool).stat().st_mode & 0o111)
]
if non_executable_bin_tools:
    raise SystemExit(
        "FAIL: modules directly invoke non-executable repository bin tools: "
        + ", ".join(non_executable_bin_tools)
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
    "reports/multiqc/multiqc_report.html",
    "reports/qc_dashboard/qc_dashboard.html",
    "qc_summary.tsv",
    "qc_summary.json",
    "top_motifs.tsv",
    "tss_profiles.tsv",
    "Dashboard data-product schemas",
    "schema version `1`",
    "generator version `1.1.0`",
    "six sequencing panels",
    "five peak panels",
    "shared 250-bp bins",
    "fragments-per-peak coverage ECDF",
    "zero-fragment peaks",
    "target palette",
    "direct endpoint labels",
    "complete AME output",
    "15 non-cognate motifs",
    "three significant digits",
    "full precision",
    "`sample_id`, `position_bp`, `signal`",
    "TSV missing numeric values are empty fields",
    "`QC` subworkflow now requires eleven inputs",
    "`TSS_ENRICHMENT.out.profiles` tuple now contains seven values",
    "qc/peaks/<sample_id>/<sample_id>.peak_qc.tsv",
    "position 0 divided by the mean first/last 100 bp",
    "visually separated",
    "no expected-motif result",
    "NA warnings rather than zero",
    "applies no biological thresholds",
]
missing_phrases = [phrase for phrase in required_phrases if phrase not in readme]
if missing_phrases:
    raise SystemExit(
        "FAIL: README is missing required documentation topics: "
        + ", ".join(missing_phrases)
    )

print(
    f"README documents all {len(config_params)} CLI parameters, "
    f"{len(environment_files)} Conda environments, manifest fields, "
    "runtime environment variables, and required operational/QC topics; "
    f"all {len(invoked_bin_tools)} module-invoked repository bin tools are executable."
)
PY
