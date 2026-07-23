#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"

required_files=(
  modules/local/validate_manifest.nf
  modules/local/demultiplex_i2.nf
  modules/local/fastqc.nf
  subworkflows/local/demultiplex.nf
)

for required_file in "${required_files[@]}"; do
  if [[ ! -f "$required_file" ]]; then
    printf 'FAIL: missing %s\n' "$required_file" >&2
    exit 1
  fi
done

for executable_script in bin/manifest.py bin/demultiplex_i2.py; do
  if [[ ! -x "$executable_script" ]]; then
    printf 'FAIL: Nextflow bin entry point is not executable: %s\n' "$executable_script" >&2
    exit 1
  fi
done

python3 - <<'PY'
from pathlib import Path

root = Path.cwd()
validate = (root / "modules/local/validate_manifest.nf").read_text()
demux = (root / "modules/local/demultiplex_i2.nf").read_text()
fastqc = (root / "modules/local/fastqc.nf").read_text()
subworkflow = (root / "subworkflows/local/demultiplex.nf").read_text()

checks = {
    "VALIDATE_MANIFEST process is declared": "process VALIDATE_MANIFEST" in validate,
    "manifest validator uses the project-qualified Python environment":
        'conda "${projectDir}/envs/python.yml"' in validate,
    "manifest validator has a pinned Python container":
        "container 'python:3.12.3-slim-bookworm'" in validate,
    "manifest validator invokes the existing CLI":
        "manifest.py validate" in validate and "${projectDir}/bin/" not in validate,
    "manifest validator stages any source basename under a fixed safe alias":
        "path manifest, stageAs: 'input_manifest.csv'" in validate,
    "manifest validator shell uses only the fixed staged manifest alias":
        '--input "input_manifest.csv"' in validate
        and "${manifest}" not in validate.split("script:", 1)[1],
    "manifest validator emits normalized JSON and versions":
        "normalized_manifest.json" in validate and "versions" in validate,
    "DEMULTIPLEX_I2 process is declared": "process DEMULTIPLEX_I2" in demux,
    "splitter uses the project-qualified Python environment":
        'conda "${projectDir}/envs/python.yml"' in demux,
    "splitter has a pinned Python container":
        "container 'python:3.12.3-slim-bookworm'" in demux,
    "splitter invokes the existing CLI":
        "demultiplex_i2.py" in demux and "${projectDir}/bin/" not in demux,
    "splitter forwards repeated samples and mismatch policy":
        "--sample" in demux and "--max-mismatches" in demux and "--allow-empty" in demux,
    "splitter validates and coerces mismatch count before interpolation":
        "Integer.parseInt" in demux
        and demux.index("Integer.parseInt") < demux.index("--max-mismatches"),
    "splitter stages synchronized inputs under collision-proof aliases":
        all(alias in demux for alias in (
            "path(r1, stageAs: 'input_R1.fastq.gz')",
            "path(r2, stageAs: 'input_R2.fastq.gz')",
            "path(i2, stageAs: 'input_I2.fastq.gz')",
        )),
    "splitter uses fixed staged and metrics filenames in the shell command":
        all(token in demux for token in (
            '--r1 "input_R1.fastq.gz"',
            '--r2 "input_R2.fastq.gz"',
            '--i2 "input_I2.fastq.gz"',
            '--metrics "demultiplex.metrics.json"',
        ))
        and '${meta.library_id}.demultiplex.metrics' not in demux,
    "splitter emits paired derived FASTQs, metrics, and versions":
        "*_R1.fastq.gz" in demux and "*_R2.fastq.gz" in demux
        and ".metrics.json" in demux and ".metrics.tsv" in demux and "versions" in demux,
    "FASTQC process is declared": "process FASTQC" in fastqc,
    "FastQC uses the project-qualified FastQC environment":
        'conda "${projectDir}/envs/fastqc.yml"' in fastqc,
    "FastQC has the pinned BioContainer":
        "container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'" in fastqc,
    "FastQC emits HTML, ZIP, and versions":
        "*_fastqc.html" in fastqc and "*_fastqc.zip" in fastqc and "versions" in fastqc,
    "subworkflow includes all Task 6 modules":
        all(name in subworkflow for name in ("VALIDATE_MANIFEST", "DEMULTIPLEX_I2", "FASTQC")),
    "subworkflow parses normalized JSON":
        "JsonSlurper" in subworkflow,
    "subworkflow groups records by physical library before demultiplexing":
        "groupTuple" in subworkflow
        and subworkflow.index("groupTuple") < subworkflow.index("DEMULTIPLEX_I2("),
    "subworkflow expands library outputs into derived sample tuples":
        "flatMap" in subworkflow and "tuple(meta, r1, r2)" in subworkflow,
    "subworkflow runs FastQC on derived sample tuples":
        "FASTQC(derived_reads)" in subworkflow,
    "subworkflow exposes reads, metrics, FastQC, and versions":
        all(token in subworkflow for token in
            ("reads = derived_reads", "metrics =", "fastqc =", "versions =")),
}

failed = [description for description, passed in checks.items() if not passed]
if failed:
    raise SystemExit("FAIL: structural assertions failed:\n  - " + "\n  - ".join(failed))
PY

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/nanocut-demultiplex.XXXXXX")
trap 'rm -rf -- "${tmp_dir:?}"' EXIT

write_fastq() {
  local output=$1
  local mate=$2
  cat > "$output" <<EOF
@read-a/${mate}
ACGT
+
IIII
@read-b/${mate}
TGCA
+
IIII
EOF
}

mkdir -p "$tmp_dir/r1-source" "$tmp_dir/r2-source" "$tmp_dir/i2-source"
write_fastq "$tmp_dir/r1-source/reads.fastq" 1
write_fastq "$tmp_dir/r2-source/reads.fastq" 2
cat > "$tmp_dir/i2-source/reads.fastq" <<'EOF'
@read-a/3
TATAGCCT
+
IIIIIIII
@read-b/3
ATAGAGGC
+
IIIIIIII
EOF

cat > "$tmp_dir/samples.csv" <<EOF
sample_id,library_id,input_group,barcode,assay_target,is_control,control_id,expected_motif,r1,r2,i2
lib_IgG,lib,25K,TATAGCCT,IgG,true,,,r1-source/reads.fastq,r2-source/reads.fastq,i2-source/reads.fastq
lib_CTCF,lib,25K,ATAGAGGC,CTCF,false,lib_IgG,CTCF,r1-source/reads.fastq,r2-source/reads.fastq,i2-source/reads.fastq
EOF

python3 bin/manifest.py validate \
  --input "$tmp_dir/samples.csv" \
  --output "$tmp_dir/normalized.json"

# A source basename containing command-substitution syntax must disappear at
# the staging boundary.  The module's fixed stageAs alias and literal shell
# argument are asserted above; this fixture exercises the corresponding file
# semantics without ever evaluating the hostile basename.
hostile_manifest="$tmp_dir/\$(touch hostile-basename-executed).csv"
safe_stage_dir="$tmp_dir/safe-manifest-stage"
cp "$tmp_dir/samples.csv" "$hostile_manifest"
mkdir -p "$safe_stage_dir"
ln -s "$hostile_manifest" "$safe_stage_dir/input_manifest.csv"
(
  cd "$safe_stage_dir"
  "$repo_root/bin/manifest.py" validate \
    --input "input_manifest.csv" \
    --output normalized_manifest.json
)
if [[ -e "$safe_stage_dir/hostile-basename-executed" ]]; then
  printf 'FAIL: hostile manifest source basename was evaluated as shell syntax\n' >&2
  exit 1
fi

mkdir -p "$tmp_dir/direct-demux"
python3 bin/demultiplex_i2.py \
  --r1 "$tmp_dir/r1-source/reads.fastq" \
  --r2 "$tmp_dir/r2-source/reads.fastq" \
  --i2 "$tmp_dir/i2-source/reads.fastq" \
  --sample lib_IgG=TATAGCCT \
  --sample lib_CTCF=ATAGAGGC \
  --max-mismatches 0 \
  --outdir "$tmp_dir/direct-demux" \
  --metrics "$tmp_dir/direct-demux/demultiplex.metrics.json"

python3 - "$tmp_dir/direct-demux" <<'PY'
import gzip
import json
import sys
from pathlib import Path

outdir = Path(sys.argv[1])
metrics = json.loads((outdir / "demultiplex.metrics.json").read_text())
assert metrics["total_reads"] == 2, metrics
assert metrics["assigned_reads"] == 2, metrics
assert metrics["ambiguous_reads"] == 0, metrics
assert metrics["unassigned_reads"] == 0, metrics
assert metrics["assignment_counts"] == {"lib_CTCF": 1, "lib_IgG": 1}, metrics
for sample_id in ("lib_IgG", "lib_CTCF"):
    for read in ("R1", "R2"):
        path = outdir / f"{sample_id}_{read}.fastq.gz"
        assert path.is_file(), path
        with gzip.open(path, "rt") as handle:
            assert sum(1 for _ in handle) == 4, path
PY

if ! command -v nextflow >/dev/null 2>&1; then
  printf 'SKIP: Nextflow runtime unavailable; static DSL2 and direct synthetic demultiplexer checks passed, but no Nextflow execution was run.\n'
  exit 0
fi

fake_bin="$tmp_dir/bin"
mkdir -p "$fake_bin"
cat > "$fake_bin/fastqc" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --version ]]; then
  printf 'FastQC v0.12.1\n'
  exit 0
fi
outdir=
inputs=()
while (( $# > 0 )); do
  case "$1" in
    --outdir)
      outdir=$2
      shift 2
      ;;
    --threads)
      shift 2
      ;;
    *)
      inputs+=("$1")
      shift
      ;;
  esac
done
mkdir -p "$outdir"
for input in "${inputs[@]}"; do
  name=$(basename "$input")
  name=${name%.gz}
  name=${name%.fastq}
  name=${name%.fq}
  printf '<html>synthetic FastQC</html>\n' > "$outdir/${name}_fastqc.html"
  : > "$outdir/${name}_fastqc.zip"
done
EOF
chmod +x "$fake_bin/fastqc"

harness="$repo_root/.nextflow-task6-test-$$.nf"
runtime_config="$tmp_dir/nextflow.config"
trace_file="$tmp_dir/trace.tsv"
trap 'rm -f -- "${harness:?}"; rm -rf -- "${tmp_dir:?}"' EXIT

cat > "$harness" <<'EOF'
nextflow.enable.dsl = 2

include { DEMULTIPLEX } from './subworkflows/local/demultiplex'

workflow {
    manifest_ch = Channel.fromPath(params.input, checkIfExists: true)
    DEMULTIPLEX(manifest_ch, params.barcode_mismatches, params.allow_empty)
    DEMULTIPLEX.out.reads.view { meta, r1, r2 ->
        "DERIVED:${meta.sample_id}:${r1.name}:${r2.name}"
    }
}
EOF

cat > "$runtime_config" <<EOF
params.outdir = '$tmp_dir/results'
conda.enabled = false
docker.enabled = false
report.enabled = false
timeline.enabled = false
trace.enabled = false
dag.enabled = false
process.executor = 'local'
process.maxForks = 1
EOF

PATH="$fake_bin:$PATH" nextflow run "$harness" \
  -c "$runtime_config" \
  --input "$tmp_dir/samples.csv" \
  --barcode_mismatches 0 \
  --allow_empty false \
  -with-trace "$trace_file" \
  -ansi-log false

demux_tasks=$(awk -F '\t' 'NR > 1 && $3 ~ /DEMULTIPLEX_I2/ { count++ } END { print count + 0 }' "$trace_file")
fastqc_tasks=$(awk -F '\t' 'NR > 1 && $3 ~ /FASTQC/ { count++ } END { print count + 0 }' "$trace_file")
[[ $demux_tasks -eq 1 ]] || {
  printf 'FAIL: expected one DEMULTIPLEX_I2 task for one physical library, observed %s\n' "$demux_tasks" >&2
  exit 1
}
[[ $fastqc_tasks -eq 2 ]] || {
  printf 'FAIL: expected two FASTQC tasks for two derived samples, observed %s\n' "$fastqc_tasks" >&2
  exit 1
}

python3 - "$tmp_dir/results" <<'PY'
import json
import sys
from pathlib import Path

results = Path(sys.argv[1])
metrics = json.loads(
    (results / "demultiplex/lib/demultiplex.metrics.json").read_text()
)
assert metrics["assignment_counts"] == {"lib_CTCF": 1, "lib_IgG": 1}, metrics
for sample_id in ("lib_IgG", "lib_CTCF"):
    for read in ("R1", "R2"):
        assert (
            results / f"demultiplex/lib/{sample_id}_{read}.fastq.gz"
        ).is_file()
        assert (
            results / f"fastqc/{sample_id}/{sample_id}_{read}_fastqc.html"
        ).is_file()
PY

printf 'PASS: one physical-library demultiplex task produced both derived samples and two per-sample FastQC tasks.\n'
