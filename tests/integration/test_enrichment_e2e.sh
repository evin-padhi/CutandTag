#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"

base_python=$(command -v python3 || true)
if [[ -z "$base_python" ]]; then
  printf 'FAIL: python3 is required for enrichment fixture setup\n' >&2
  exit 1
fi

plot_python=""
for candidate in python3.11 python3.12 python3 python; do
  if ! command -v "$candidate" >/dev/null 2>&1; then
    continue
  fi
  if "$candidate" - <<'PY' >/dev/null 2>&1
import matplotlib
PY
  then
    plot_python=$(command -v "$candidate")
    break
  fi
done

if [[ -z "$plot_python" ]]; then
  printf 'FAIL: no Python interpreter with matplotlib is available for real enrichment plot generation\n' >&2
  exit 1
fi

required_files=(
  main.nf
  conf/test.config
  bin/peak_enrichment.py
  tests/data/enrichment/chipseq.csv
  tests/data/enrichment/ctcf.bed
  tests/data/enrichment/runx1.bed
  tests/data/enrichment/foregrounds.tsv
  tests/data/enrichment/reference.fa
)

for required_file in "${required_files[@]}"; do
  if [[ ! -f "$required_file" ]]; then
    printf 'FAIL: missing %s\n' "$required_file" >&2
    exit 1
  fi
done

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/nanocut-enrichment-e2e.XXXXXX")
trap 'rm -rf -- "${tmp_dir:?}"' EXIT

combined_manifest="$tmp_dir/reference_manifest.tsv"
"$base_python" - "$repo_root" "$combined_manifest" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
output = Path(sys.argv[2])
foreground_manifest = root / "tests/data/enrichment/foregrounds.tsv"
chipseq_manifest = root / "tests/data/enrichment/chipseq.csv"
fixture_dir = root / "tests/data/enrichment"

with foreground_manifest.open(newline="", encoding="utf-8") as handle:
    foreground_rows = list(csv.DictReader(handle, delimiter="\t"))

with chipseq_manifest.open(newline="", encoding="utf-8") as handle:
    chipseq_rows = list(csv.DictReader(handle))

with output.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["reference_id", "tf", "reference_type", "peak_file"],
        delimiter="\t",
    )
    writer.writeheader()
    for row in chipseq_rows:
        writer.writerow(
            {
                "reference_id": row["reference_id"],
                "tf": row["tf"],
                "reference_type": "chipseq",
                "peak_file": str((fixture_dir / row["peak_file"]).resolve()),
            }
        )
    for row in foreground_rows:
        writer.writerow(
            {
                "reference_id": row["foreground_id"],
                "tf": row["foreground_tf"],
                "reference_type": "called_tf",
                "peak_file": str((fixture_dir / row["peak_file"]).resolve()),
            }
        )
PY

cli_outdir="$tmp_dir/cli"
mpl_config_dir="$tmp_dir/mplconfig"
mkdir -p "$mpl_config_dir"

MPLCONFIGDIR="$mpl_config_dir" "$plot_python" bin/peak_enrichment.py \
  --foreground-manifest tests/data/enrichment/foregrounds.tsv \
  --reference-manifest "$combined_manifest" \
  --fasta tests/data/enrichment/reference.fa \
  --outdir "$cli_outdir" \
  --permutations 20 \
  --seed 1729 \
  --gc-tolerance 0.02

"$base_python" - "$cli_outdir" <<'PY'
import csv
import math
import sys
from pathlib import Path

outdir = Path(sys.argv[1])
models = {
    "random",
    "length_matched",
    "gc_matched",
    "length_gc_matched",
}
required_outputs = [
    outdir / "peak_enrichment.tsv",
    outdir / "enrichment_status.tsv",
    outdir / "observed_vs_null.png",
]
required_outputs.extend(outdir / f"matrix_{model}.tsv" for model in models)
required_outputs.extend(outdir / f"matrix_{model}.png" for model in models)

missing = [
    str(path)
    for path in required_outputs
    if not path.is_file() or path.stat().st_size == 0
]
if missing:
    raise SystemExit("FAIL: missing or empty CLI outputs:\n  - " + "\n  - ".join(missing))

with (outdir / "peak_enrichment.tsv").open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

if not rows:
    raise SystemExit("FAIL: peak_enrichment.tsv is empty")

if {row["background_model"] for row in rows} != models:
    raise SystemExit(
        "FAIL: unexpected background models: "
        + repr(sorted({row["background_model"] for row in rows}))
    )

numeric_fields = (
    "observed_overlap_count",
    "null_mean_overlap_count",
    "null_stddev_overlap_count",
    "enrichment_ratio",
    "empirical_p_value",
)
for row in rows:
    if row["status"] != "ok":
        raise SystemExit(f"FAIL: expected ok status rows, found {row['status']!r}")
    for field in numeric_fields:
        if not math.isfinite(float(row[field])):
            raise SystemExit(f"FAIL: non-finite {field} in row {row}")

with (outdir / "enrichment_status.tsv").open(newline="", encoding="utf-8") as handle:
    status_rows = list(csv.DictReader(handle, delimiter="\t"))

if len(status_rows) != len(rows):
    raise SystemExit(
        f"FAIL: expected {len(rows)} status rows, found {len(status_rows)}"
    )

if {row["background_model"] for row in status_rows} != models:
    raise SystemExit("FAIL: status table background models did not match the enrichment rows")

png_signature = b"\x89PNG\r\n\x1a\n"
for path in sorted(outdir.glob("*.png")):
    data = path.read_bytes()
    if not data.startswith(png_signature):
        raise SystemExit(f"FAIL: {path.name} is not a PNG file")
    if len(data) <= 500:
        raise SystemExit(f"FAIL: {path.name} is too small to be a real matplotlib plot")
PY

malformed_manifest="$tmp_dir/malformed.csv"
"$base_python" - "$malformed_manifest" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
with path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["reference_id", "tf", "peak_file"])
    writer.writeheader()
    fixture_dir = Path.cwd() / "tests/data/enrichment"
    writer.writerow({"reference_id": "dup", "tf": "CTCF", "peak_file": str((fixture_dir / "ctcf.bed").resolve())})
    writer.writerow({"reference_id": "dup", "tf": "RUNX1", "peak_file": str((fixture_dir / "runx1.bed").resolve())})
PY

set +e
"$base_python" bin/peak_enrichment.py \
  --foreground-manifest tests/data/enrichment/foregrounds.tsv \
  --reference-manifest "$malformed_manifest" \
  --fasta tests/data/enrichment/reference.fa \
  --outdir "$tmp_dir/malformed-out" \
  --permutations 20 \
  --seed 1729 \
  --gc-tolerance 0.02 \
  >"$tmp_dir/malformed.stdout" 2>"$tmp_dir/malformed.stderr"
status=$?
set -e

if [[ $status -eq 0 ]]; then
  printf 'FAIL: malformed manifest unexpectedly passed the CLI\n' >&2
  exit 1
fi

if ! grep -q 'duplicate reference_id' "$tmp_dir/malformed.stderr"; then
  printf 'FAIL: malformed manifest did not report duplicate reference_id\n' >&2
  exit 1
fi

if ! "$base_python" - <<'PY' >/dev/null 2>&1
import matplotlib
PY
then
  set +e
  "$base_python" bin/peak_enrichment.py \
    --foreground-manifest tests/data/enrichment/foregrounds.tsv \
    --reference-manifest "$combined_manifest" \
    --fasta tests/data/enrichment/reference.fa \
    --outdir "$tmp_dir/base-python-out" \
    --permutations 20 \
    --seed 1729 \
    --gc-tolerance 0.02 \
    >"$tmp_dir/base-python.stdout" 2>"$tmp_dir/base-python.stderr"
  missing_dep_status=$?
  set -e

  if [[ $missing_dep_status -eq 0 ]]; then
    printf 'FAIL: matplotlib-missing interpreter unexpectedly generated plots successfully\n' >&2
    exit 1
  fi

  if ! grep -qi 'matplotlib' "$tmp_dir/base-python.stderr"; then
    printf 'FAIL: matplotlib-missing interpreter did not report the missing plotting dependency\n' >&2
    exit 1
  fi
fi

if ! command -v nextflow >/dev/null 2>&1; then
  printf '%s\n' \
    'SKIP: direct offline enrichment CLI fixture checks passed, but Nextflow runtime is unavailable so the pipeline/MultiQC enrichment run was not exercised.'
  exit 0
fi

runtime_dir="$tmp_dir/runtime"
mkdir -p "$runtime_dir"

(
  cd "$runtime_dir"
  nextflow run "$repo_root/main.nf" \
    -profile test \
    --outdir "$runtime_dir/results" \
    -work-dir "$runtime_dir/work" \
    --fasta "$repo_root/tests/data/enrichment/reference.fa" \
    --chipseq_input "$repo_root/tests/data/enrichment/chipseq.csv" \
    --enrichment_permutations 20 \
    --enrichment_seed 1729 \
    --enrichment_gc_tolerance 0.02
) >"$tmp_dir/nextflow.log" 2>&1

python3 - "$runtime_dir/results" <<'PY'
import csv
import sys
from pathlib import Path

results = Path(sys.argv[1])
required_outputs = [
    results / "peak_enrichment/peak_enrichment.tsv",
    results / "peak_enrichment/enrichment_status.tsv",
    results / "peak_enrichment/plots/observed_vs_null.png",
    results / "reports/summary/peak_enrichment.tsv",
    results / "reports/multiqc/multiqc_report.html",
    results / "reports/multiqc/observed_vs_null.png",
]
required_outputs.extend(
    results / "peak_enrichment/plots" / f"matrix_{model}.tsv"
    for model in ("random", "length_matched", "gc_matched", "length_gc_matched")
)
required_outputs.extend(
    results / "peak_enrichment/plots" / f"matrix_{model}.png"
    for model in ("random", "length_matched", "gc_matched", "length_gc_matched")
)
required_outputs.extend(
    results / "reports/summary" / f"matrix_{model}.tsv"
    for model in ("random", "length_matched", "gc_matched", "length_gc_matched")
)
required_outputs.extend(
    results / "reports/multiqc" / f"matrix_{model}.png"
    for model in ("random", "length_matched", "gc_matched", "length_gc_matched")
)

missing = [
    str(path)
    for path in required_outputs
    if not path.is_file() or path.stat().st_size == 0
]
if missing:
    raise SystemExit(
        "FAIL: missing or empty enrichment pipeline outputs:\n  - "
        + "\n  - ".join(missing)
    )

run_summary = {
    row["metric"]: row["value"]
    for row in csv.DictReader(
        (results / "pipeline_info/run_summary.txt").open(
            newline="", encoding="utf-8"
        ),
        delimiter="\t",
        fieldnames=["metric", "value"],
    )
}
if run_summary.get("peak_enrichment") != "reports/summary/peak_enrichment.tsv":
    raise SystemExit(
        "FAIL: run summary did not advertise the enrichment summary output"
    )
PY
