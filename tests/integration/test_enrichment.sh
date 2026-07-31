#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"

required_files=(
  main.nf
  nextflow.config
  nextflow_schema.json
  conf/test.config
  modules/local/validate_chipseq_manifest.nf
  subworkflows/local/enrichment.nf
)

for required_file in "${required_files[@]}"; do
  if [[ ! -f "$required_file" ]]; then
    printf 'FAIL: missing %s\n' "$required_file" >&2
    exit 1
  fi
done

python3 - <<'PY'
import json
import re
from pathlib import Path

root = Path.cwd()
main = (root / "main.nf").read_text(encoding="utf-8")
base_config = (root / "nextflow.config").read_text(encoding="utf-8")
test_config = (root / "conf/test.config").read_text(encoding="utf-8")
schema = json.loads((root / "nextflow_schema.json").read_text(encoding="utf-8"))
validate_chipseq = (root / "modules/local/validate_chipseq_manifest.nf").read_text(
    encoding="utf-8"
)
enrichment = (root / "subworkflows/local/enrichment.nf").read_text(
    encoding="utf-8"
)

props = schema["properties"]
checks = {
    "main includes the ENRICHMENT subworkflow":
        re.search(r"include\s+\{\s*ENRICHMENT\s*\}\s+from\s+'\.\/subworkflows/local/enrichment'", main),
    "chipseq_input is validated as an optional readable CSV":
        "--chipseq_input" in main
        and "validated.chipseq_input" in main
        and "validateRegularFile(" in main
        and ".endsWith('.csv')" in main,
    "enrichment parameters are validated and surfaced in validated_parameters.json":
        all(token in main for token in (
            "validated.enrichment_permutations",
            "validated.enrichment_seed",
            "validated.enrichment_gc_tolerance",
            "--enrichment_permutations",
            "--enrichment_seed",
            "--enrichment_gc_tolerance",
        ))
        and "parameter_map = new LinkedHashMap(validated)" in main,
    "peak enrichment requires FASTA only when chipseq_input is supplied":
        "chipseq_input" in main
        and re.search(
            r"validated\.chipseq_input\s*!=\s*null\s*&&\s*validated\.fasta\s*==\s*null",
            main,
        )
        and "--chipseq_input is supplied" in main,
    "new runtime defaults are defined in the base config":
        all(token in base_config for token in (
            "chipseq_input            = null",
            "enrichment_permutations  = 1000",
            "enrichment_seed          = 1729",
            "enrichment_gc_tolerance  = 0.02",
        )),
    "test profile leaves chipseq_input disabled":
        "chipseq_input = null" in test_config,
    "schema exposes the optional enrichment parameters and defaults":
        props["chipseq_input"]["type"] == "string"
        and props["chipseq_input"]["format"] == "file-path"
        and props["enrichment_permutations"]["default"] == 1000
        and props["enrichment_seed"]["default"] == 1729
        and props["enrichment_gc_tolerance"]["default"] == 0.02,
    "existing required-input rules remain intact":
        schema["required"] == ["input", "macs_genome_size"]
        and schema["allOf"],
    "main constructs an empty chipseq channel when the manifest is absent":
        "chipseq_manifest_ch = optionalPathChannel(validated.chipseq_input)" in main
        or "chipseq_manifest_ch = Channel.empty()" in main,
    "main wires ENRICHMENT conditionally and merges its versions":
        "ENRICHMENT(" in main
        and "enrichment_versions_ch = Channel.empty()" in main
        and "enrichment_results_ch = Channel.empty()" in main
        and "enrichment_status_ch = Channel.empty()" in main
        and "enrichment_plots_ch = Channel.empty()" in main
        and "ENRICHMENT.out.versions" in main,
    "chipseq manifest validation module is declared":
        "process VALIDATE_CHIPSEQ_MANIFEST" in validate_chipseq
        and 'conda "${projectDir}/envs/python.yml"' in validate_chipseq
        and "peak_enrichment.py" in validate_chipseq
        and "normalized_chipseq_manifest.tsv" in validate_chipseq,
    "ENRICHMENT validates metadata and prepares normalized manifests":
        "workflow ENRICHMENT" in enrichment
        and "foreground_id" in enrichment
        and "reference_id" in enrichment
        and "reference_type" in enrichment
        and "called_tf" in enrichment
        and "VALIDATE_CHIPSEQ_MANIFEST" in enrichment,
    "ENRICHMENT publishes results, status, plots, and versions":
        all(token in enrichment for token in (
            "emit:",
            "results =",
            "plots =",
            "status =",
            "versions =",
            "peak_enrichment.tsv",
            "enrichment_status.tsv",
            "observed_vs_null.png",
        )),
}

failed = [description for description, passed in checks.items() if not passed]
if failed:
    raise SystemExit(
        "FAIL: enrichment structural assertions failed:\n  - "
        + "\n  - ".join(failed)
    )
PY
