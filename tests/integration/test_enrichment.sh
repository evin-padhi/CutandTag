#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"

required_files=(
  main.nf
  nextflow.config
  nextflow_schema.json
  conf/test.config
  bin/peak_enrichment.py
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
peak_enrichment = (root / "bin/peak_enrichment.py").read_text(encoding="utf-8")
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
    "launch-time chipseq validation uses the tested Python CSV utility":
        "validateChipseqManifestAtLaunch" in main
        and "validate-chipseq-manifest" in main
        and "new ProcessBuilder" in main
        and re.search(
            r"validated\.chipseq_input\s*=\s*validateChipseqManifestAtLaunch\(",
            main,
        ),
    "peak enrichment requires FASTA only when chipseq_input is supplied":
        "chipseq_input" in main
        and "peak enrichment requires --fasta when --chipseq_input is supplied" in main
        and "validateChipseqManifestAtLaunch" in main,
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
    "main parses public chipseq CSV with splitCsv and forces public references to chipseq":
        ".splitCsv(header: true)" in main
        and "reference_type: 'chipseq'" in main
        and "chipseq_reference_inputs_ch = Channel.empty()" in main
        and "['reference_type']" not in main,
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
        and "normalized_chipseq_manifest.tsv" in validate_chipseq
        and "csv.DictWriter" in validate_chipseq,
    "chipseq manifest validation tag does not reference the removed manifest input":
        "tag 'chipseq-manifest'" in validate_chipseq
        and "manifest.simpleName" not in validate_chipseq,
    "launch and staged manifest validation use the streaming FASTA size loader":
        "load_fasta_chrom_sizes(args.fasta)" in peak_enrichment
        and "load_fasta_chrom_sizes" in validate_chipseq
        and "load_fasta_sequences" not in validate_chipseq,
    "every external chipseq BED is a staged path dependency and emitted for enrichment":
        "path(reference_peak_files" in validate_chipseq
        and 'path "reference_peaks/*.bed*", emit: peaks' in validate_chipseq
        and enrichment.count("path external_peak_files") >= 2
        and "VALIDATE_CHIPSEQ_MANIFEST.out.peaks" in enrichment,
    "gzip chipseq references retain their suffix after staging":
        'f"{record.reference_id}.bed.gz"' in validate_chipseq
        and 'record.peak_file.name.lower().endswith(".gz")' in validate_chipseq
        and 'path "reference_peaks/*.bed*", emit: peaks' in validate_chipseq,
    "normalized chipseq manifests contain only staged relative peak paths":
        '"peak_file": f"reference_peaks/{destination.name}"' in validate_chipseq
        and "str(record.peak_file)" not in validate_chipseq
        and "chipseq/reference_manifest.tsv" in enrichment,
    "repository enrichment scripts are staged safely instead of source-interpolated":
        "path enrichment_script, stageAs: 'bin/peak_enrichment.py'" in validate_chipseq
        and enrichment.count("path enrichment_script, stageAs: 'bin/peak_enrichment.py'") >= 2
        and "projectBin" not in validate_chipseq
        and "projectBin" not in enrichment,
    "Docker enrichment uses a pinned Python image containing matplotlib":
        all(
            re.search(
                rf"process {process_name} \{{.*?"
                r"container 'quay.io/jupyter/scipy-notebook:82d322f00937'",
                enrichment,
                re.S,
            )
            for process_name in (
                "RUN_PEAK_ENRICHMENT",
                "WRITE_EMPTY_ENRICHMENT_OUTPUTS",
            )
        ),
    "called-TF identifiers are collision-resistant and DSL2 wiring avoids deprecated into":
        "MessageDigest.getInstance('SHA-256')" in enrichment
        and "foregroundHash" in enrichment
        and ".into {" not in enrichment,
    "enrichment artifacts publish under the specified enrichment directory":
        enrichment.count('${params.outdir}/enrichment') >= 2
        and '${params.outdir}/peak_enrichment' not in enrichment,
    "ENRICHMENT validates metadata and prepares normalized manifests":
        "workflow ENRICHMENT" in enrichment
        and "foreground_id" in enrichment
        and "reference_id" in enrichment
        and "reference_type" in enrichment
        and "called_tf" in enrichment
        and "VALIDATE_CHIPSEQ_MANIFEST" in enrichment,
    "ENRICHMENT has an explicit zero-foreground status-only branch":
        "WRITE_EMPTY_ENRICHMENT_OUTPUTS" in enrichment
        and "rows.isEmpty()" in enrichment
        and 'foreground_manifest.tsv' in enrichment,
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

if ! command -v nextflow >/dev/null 2>&1; then
  printf '%s\n' \
    'SKIP: Nextflow runtime unavailable; static enrichment wiring checks passed, but launch-time malformed-manifest execution was not exercised.'
  exit 0
fi

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/nanocut-enrichment.XXXXXX")
trap 'rm -rf -- "${tmp_dir:?}"' EXIT

cat > "$tmp_dir/duplicate.bed" <<'EOF'
chrMini	10	20
EOF

cat > "$tmp_dir/chipseq.csv" <<EOF
reference_id,tf,peak_file
dup,CTCF,$tmp_dir/duplicate.bed
dup,RUNX1,$tmp_dir/duplicate.bed
EOF

set +e
nextflow run main.nf -profile test -stub-run --chipseq_input "$tmp_dir/chipseq.csv" \
  >"$tmp_dir/launch-validation.log" 2>&1
status=$?
set -e

if [[ $status -eq 0 ]]; then
  printf 'FAIL: malformed chipseq manifest unexpectedly passed launch-time validation\n' >&2
  exit 1
fi

if ! grep -q 'duplicate reference_id' "$tmp_dir/launch-validation.log"; then
  printf 'FAIL: malformed chipseq manifest did not report duplicate reference_id at launch time\n' >&2
  exit 1
fi

if grep -q 'MACS2_BROAD\|DEMULTIPLEX\|executor >' "$tmp_dir/launch-validation.log"; then
  printf 'FAIL: malformed chipseq manifest reached workflow scheduling before validation stopped the run\n' >&2
  exit 1
fi
