#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"

required_files=(
  main.nf
  conf/test.config
  tests/data/e2e/samples.csv
  tests/data/e2e/reads_R1.fastq
  tests/data/e2e/reads_R2.fastq
  tests/data/e2e/reads_I2.fastq
  tests/data/e2e/reference.fa
  tests/data/e2e/blacklist.bed
  tests/data/e2e/motifs.meme
  tests/data/e2e/fakebin/fake_bio_tool.py
)

for required_file in "${required_files[@]}"; do
  if [[ ! -f "$required_file" ]]; then
    printf 'FAIL: missing %s\n' "$required_file" >&2
    exit 1
  fi
done

python3 - <<'PY'
import csv
import re
from pathlib import Path

root = Path.cwd()
main = (root / "main.nf").read_text(encoding="utf-8")
config = (root / "conf/test.config").read_text(encoding="utf-8")
base_config = (root / "nextflow.config").read_text(encoding="utf-8")

checks = {
    "main includes all five subworkflows":
        all(
            re.search(rf"include\s+\{{\s*{name}\s*\}}", main)
            for name in ("DEMULTIPLEX", "ALIGN_QC", "PEAKS", "MOTIFS", "QC")
        ),
    "workflow calls follow complete dependency order":
        all(
            main.index(first) < main.index(second)
            for first, second in (
                ("\n    DEMULTIPLEX(", "\n    ALIGN_QC("),
                ("\n    ALIGN_QC(", "\n    PEAKS("),
                ("\n    PEAKS(", "\n        MOTIFS("),
                ("\n        MOTIFS(", "\n    QC("),
            )
        ),
    "demultiplex call matches its three-input contract":
        re.search(
            r"DEMULTIPLEX\(\s*manifest_ch,\s*barcode_mismatches_ch,\s*allow_empty_ch\s*\)",
            main,
            re.S,
        ),
    "alignment call matches its four-input contract":
        re.search(
            r"ALIGN_QC\(\s*DEMULTIPLEX\.out\.reads,\s*validated\.fasta,\s*"
            r"validated\.bowtie2_index,\s*min_mapq_ch\s*\)",
            main,
            re.S,
        ),
    "peak call matches its four-input contract":
        re.search(
            r"PEAKS\(\s*ALIGN_QC\.out\.analysis_bam,\s*blacklist_ch,\s*"
            r"macs_genome_size_ch,\s*narrow_peaks_ch\s*\)",
            main,
            re.S,
        ),
    "motif call matches its exact seven-input contract in declared order":
        re.search(
            r"MOTIFS\(\s*PEAKS\.out\.final_broad_peaks,\s*"
            r"PEAKS\.out\.motif_summits,\s*motif_fasta_ch,\s*"
            r"motif_blacklist_ch,\s*motif_db_ch,\s*motif_window_ch,\s*"
            r"motif_use_narrow_peaks_ch\s*\)",
            main,
            re.S,
        ),
    "QC call matches its nine-input contract":
        re.search(
            r"QC\(\s*ALIGN_QC\.out\.filtered_bam,\s*"
            r"PEAKS\.out\.final_broad_peaks,\s*ALIGN_QC\.out\.coverage,\s*"
            r"ALIGN_QC\.out\.metrics,\s*DEMULTIPLEX\.out\.metrics,\s*"
            r"DEMULTIPLEX\.out\.fastqc,\s*motif_metrics_ch,\s*gtf_ch,\s*tss_bed_ch\s*\)",
            main,
            re.S,
        ),
    "disabled motif analysis produces closed empty channels":
        "motif_metrics_ch = Channel.empty()" in main
        and "motif_versions_ch = Channel.empty()" in main,
    "optional paths are represented by closed channels":
        "optionalPathChannel" in main and "Channel.empty()" in main,
    "manifest and regular-file paths are validated before channels are created":
        "validateRegularFile" in main
        and "validated.input" in main
        and main.index("validatePipelineParameters") < main.index("Channel.fromPath"),
    "reference validation requires FASTA or a complete Bowtie2 index":
        "requires --fasta or --bowtie2_index" in main
        and "validateBowtie2IndexPrefix" in main
        and all(suffix in main for suffix in ("rev.1", "rev.2", ".bt2l")),
    "motif analysis requires FASTA and validates MEME syntax":
        "requires --fasta when --motif_db is supplied" in main
        and "validateMemeDatabase" in main
        and "MEME version" in main,
    "fixed NanoScope MACS settings are workflow-validated":
        all(name in main for name in (
            "macs_llocal",
            "macs_keep_dup",
            "macs_broad_cutoff",
            "macs_max_gap",
        ))
        and "must remain" in main,
    "barcode, MAPQ, boolean, and motif settings are strongly validated":
        all(token in main for token in (
            "validateNonNegativeInteger",
            "validateMapq",
            "validateBoolean",
            "validatePositiveInteger",
            "barcode_mismatches",
            "allow_empty",
            "motif_use_narrow_peaks",
            "motif_window",
        )),
    "all subworkflow versions are merged and collected":
        all(token in main for token in (
            "DEMULTIPLEX.out.versions",
            "ALIGN_QC.out.versions",
            "PEAKS.out.versions",
            "motif_versions_ch",
            "QC.out.versions",
            "COLLECT_VERSIONS",
            "software_versions.yml",
        )),
    "all major metrics are merged and emitted":
        all(token in main for token in (
            "DEMULTIPLEX.out.metrics",
            "ALIGN_QC.out.metrics",
            "QC.out.target_qc",
            "motif_metrics_ch",
            "metrics = all_metrics_ch",
        )),
    "pipeline info retains validated parameters and a completion summary":
        all(token in main for token in (
            "validated_parameters.json",
            "run_summary.txt",
            "WRITE_PIPELINE_PARAMETERS",
            "WRITE_COMPLETION_SUMMARY",
            "workflow.onComplete",
        )),
    "test profile is included from the base configuration":
        re.search(
            r"test\s*\{[^}]*includeConfig\s+['\"]conf/test\.config['\"]",
            base_config,
            re.S,
        ),
    "test profile is local, serial, and disables dependency/network backends":
        all(token in config for token in (
            "executor = 'local'",
            "maxForks = 1",
            "conda.enabled = false",
            "docker.enabled = false",
            "singularity.enabled = false",
        )),
    "test profile uses only repository fixtures and fake tools":
        all(token in config for token in (
            "tests/data/e2e/samples.csv",
            "tests/data/e2e/reference.fa",
            "tests/data/e2e/blacklist.bed",
            "tests/data/e2e/motifs.meme",
            "tests/data/e2e/fakebin",
        )),
}

failed = [description for description, passed in checks.items() if not passed]
if failed:
    raise SystemExit(
        "FAIL: end-to-end structural assertions failed:\n  - "
        + "\n  - ".join(failed)
    )

with (root / "tests/data/e2e/samples.csv").open(
    newline="", encoding="utf-8"
) as handle:
    rows = list(csv.DictReader(handle))
assert len(rows) == 2
assert {row["barcode"] for row in rows} == {"TATAGCCT", "ATAGAGGC"}
controls = [row for row in rows if row["is_control"].lower() == "true"]
targets = [row for row in rows if row["is_control"].lower() == "false"]
assert len(controls) == len(targets) == 1
assert targets[0]["control_id"] == controls[0]["sample_id"]
assert targets[0]["expected_motif"] == "CTCF"
assert targets[0]["input_group"] == controls[0]["input_group"]

reference = (root / "tests/data/e2e/reference.fa").read_text(encoding="utf-8")
assert reference.startswith(">chrMini\n")
assert len("".join(reference.splitlines()[1:])) >= 4000
motifs = (root / "tests/data/e2e/motifs.meme").read_text(encoding="utf-8")
assert motifs.startswith("MEME version ")
assert "MOTIF MA0139.1 CTCF" in motifs
PY

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/nanocut-e2e.XXXXXX")
trap 'rm -rf -- "${tmp_dir:?}"' EXIT

python3 bin/manifest.py validate \
  --input tests/data/e2e/samples.csv \
  --output "$tmp_dir/normalized.json"

python3 bin/demultiplex_i2.py \
  --r1 tests/data/e2e/reads_R1.fastq \
  --r2 tests/data/e2e/reads_R2.fastq \
  --i2 tests/data/e2e/reads_I2.fastq \
  --sample MINI_IgG=TATAGCCT \
  --sample MINI_CTCF=ATAGAGGC \
  --max-mismatches 0 \
  --outdir "$tmp_dir/demux" \
  --metrics "$tmp_dir/demultiplex.metrics.json"

python3 - "$tmp_dir" <<'PY'
import gzip
import json
import sys
from pathlib import Path

tmp = Path(sys.argv[1])
with (tmp / "demultiplex.metrics.json").open(encoding="utf-8") as handle:
    metrics = json.load(handle)
assert metrics["total_reads"] == 4
assert metrics["assigned_reads"] == 4
assert metrics["ambiguous_reads"] == 0
assert metrics["unassigned_reads"] == 0
assert metrics["assignment_counts"] == {"MINI_CTCF": 2, "MINI_IgG": 2}
for sample in ("MINI_IgG", "MINI_CTCF"):
    for mate in ("R1", "R2"):
        path = tmp / "demux" / f"{sample}_{mate}.fastq.gz"
        assert path.stat().st_size > 0
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            assert sum(1 for _ in handle) == 8
PY

if ! command -v nextflow >/dev/null 2>&1; then
  printf '%s\n' \
    'SKIP: Nextflow runtime unavailable; strong static wiring/validation and direct two-barcode fixture checks passed, but no DSL2 or resume execution was run.'
  exit 0
fi

runtime_dir="$tmp_dir/runtime"
mkdir -p "$runtime_dir"

first_log="$tmp_dir/first.log"
second_log="$tmp_dir/resume.log"
first_trace="$tmp_dir/first.trace.txt"
second_trace="$tmp_dir/resume.trace.txt"

(
  cd "$runtime_dir"
  nextflow run "$repo_root/main.nf" \
    -profile test \
    --outdir "$runtime_dir/results" \
    -work-dir "$runtime_dir/work" \
    -with-trace "$first_trace"
) >"$first_log" 2>&1

python3 - "$runtime_dir/results" <<'PY'
import csv
import json
import sys
from pathlib import Path

results = Path(sys.argv[1])
required_nonempty = [
    results / "demultiplex/MINI/MINI_IgG_R1.fastq.gz",
    results / "demultiplex/MINI/MINI_CTCF_R1.fastq.gz",
    results / "alignment/MINI_IgG/analysis.bam",
    results / "alignment/MINI_CTCF/analysis.bam",
    results / "peaks/MINI_CTCF/broad/final/final.broadPeak",
    results / "qc/peaks/MINI_CTCF/MINI_CTCF.peak_qc.tsv",
    results / "motifs/MINI_CTCF/expected_motif_qc/expected_motif_qc.tsv",
    results / "reports/summary/combined_target_qc.tsv",
    results / "reports/multiqc/multiqc_report.html",
    results / "pipeline_info/software_versions.yml",
    results / "pipeline_info/validated_parameters.json",
    results / "pipeline_info/run_summary.txt",
]
for path in required_nonempty:
    assert path.is_file() and path.stat().st_size > 0, path

with (results / "demultiplex/MINI/demultiplex.metrics.json").open(
    encoding="utf-8"
) as handle:
    demux = json.load(handle)
assert demux["total_reads"] == 4
assert demux["assigned_reads"] == 4

with (results / "qc/peaks/MINI_CTCF/MINI_CTCF.peak_qc.tsv").open(
    newline="", encoding="utf-8"
) as handle:
    peak_metrics = {
        row["metric"]: row["value"]
        for row in csv.DictReader(handle, delimiter="\t")
    }
assert peak_metrics["sample_id"] == "MINI_CTCF"
assert int(peak_metrics["peak_count"]) >= 1
assert float(peak_metrics["frip"]) > 0

with (
    results
    / "motifs/MINI_CTCF/expected_motif_qc/expected_motif_qc.tsv"
).open(newline="", encoding="utf-8") as handle:
    motif_metrics = {
        row["metric"]: row["value"]
        for row in csv.DictReader(handle, delimiter="\t")
    }
assert motif_metrics["status"] == "pass"
PY

(
  cd "$runtime_dir"
  nextflow run "$repo_root/main.nf" \
    -profile test \
    --outdir "$runtime_dir/results" \
    -work-dir "$runtime_dir/work" \
    -with-trace "$second_trace" \
    -resume
) >"$second_log" 2>&1

python3 - "$first_trace" "$second_trace" <<'PY'
import csv
import sys
from collections import Counter
from pathlib import Path


def rows(path_text):
    path = Path(path_text)
    if not path.is_file() or not path.stat().st_size:
        raise SystemExit(f"FAIL: missing or empty Nextflow trace: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        result = list(csv.DictReader(handle, delimiter="\t"))
    if not result or not {"name", "status"}.issubset(result[0]):
        raise SystemExit(f"FAIL: malformed Nextflow trace: {path}")
    return result


first = rows(sys.argv[1])
resumed = rows(sys.argv[2])
bad_first = [
    (row["name"], row["status"])
    for row in first
    if row["status"].upper() != "COMPLETED"
]
if bad_first:
    raise SystemExit(
        "FAIL: first run has non-COMPLETED process rows: "
        + repr(bad_first)
    )

expected_process_types = {
    "VALIDATE_MANIFEST",
    "DEMULTIPLEX_I2",
    "FASTQC",
    "BOWTIE2_BUILD",
    "BOWTIE2_ALIGN",
    "SAMTOOLS_SORT_INDEX",
    "SAMTOOLS_FILTER",
    "SAMTOOLS_METRICS",
    "BAMCOVERAGE",
    "MACS2_BROAD",
    "MACS2_NARROW",
    "FILTER_BLACKLIST",
    "PREPARE_MOTIF_SEQUENCES",
    "AME",
    "STREME",
    "FIMO",
    "MOTIF_SUMMARY",
    "BAM_TO_FRAGMENTS",
    "FILTERED_BAM_QC",
    "PEAK_QC",
    "DEMUX_QC_CUSTOM",
    "LIBRARY_QC_CUSTOM",
    "MOTIF_QC_CUSTOM",
    "MULTIQC",
    "WRITE_PIPELINE_PARAMETERS",
    "COLLECT_VERSIONS",
    "WRITE_COMPLETION_SUMMARY",
}


def process_type(name):
    return name.rsplit(":", 1)[-1].split(" (", 1)[0]


observed_types = {process_type(row["name"]) for row in first}
missing_types = sorted(expected_process_types - observed_types)
if missing_types:
    raise SystemExit(
        "FAIL: first run omitted expected process types: "
        + ", ".join(missing_types)
    )

first_tasks = Counter(row["name"] for row in first)
resumed_cached = Counter(
    row["name"]
    for row in resumed
    if row["status"].upper() == "CACHED"
)
non_cached = [
    (row["name"], row["status"])
    for row in resumed
    if row["status"].upper() != "CACHED"
]
if non_cached:
    raise SystemExit(
        "FAIL: resumed run executed/submitted non-cached process rows: "
        + repr(non_cached)
    )
if resumed_cached != first_tasks:
    missing = first_tasks - resumed_cached
    extra = resumed_cached - first_tasks
    raise SystemExit(
        "FAIL: resumed cached task set differs from first completed task set; "
        f"missing={dict(missing)} extra={dict(extra)}"
    )
PY

printf '%s\n' \
  'PASS: synthetic demultiplex-to-report workflow executed with fake tools and every first-run process task was cached on resume with zero new task executions.'
