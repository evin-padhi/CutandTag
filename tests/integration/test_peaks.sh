#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"

required_files=(
  modules/local/macs2_broad.nf
  modules/local/macs2_narrow.nf
  modules/local/filter_blacklist.nf
  subworkflows/local/peaks.nf
)

for required_file in "${required_files[@]}"; do
  if [[ ! -f "$required_file" ]]; then
    printf 'FAIL: missing %s\n' "$required_file" >&2
    exit 1
  fi
done

python3 - <<'PY'
import csv
from pathlib import Path

root = Path.cwd()
broad = (root / "modules/local/macs2_broad.nf").read_text()
narrow = (root / "modules/local/macs2_narrow.nf").read_text()
blacklist = (root / "modules/local/filter_blacklist.nf").read_text()
peaks = (root / "subworkflows/local/peaks.nf").read_text()

macs2_container = "quay.io/biocontainers/macs2:2.2.9.1--py39hbcbf7aa_4"
bedtools_container = "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"

checks = {
    "MACS2_BROAD process is declared":
        "process MACS2_BROAD" in broad,
    "broad calling uses the pinned MACS2 environment and image":
        'conda "${projectDir}/envs/macs2.yml"' in broad
        and f"container '{macs2_container}'" in broad,
    "broad calling stages target and control BAMs separately":
        all(token in broad for token in (
            "stageAs: 'target.analysis.bam'",
            "stageAs: 'target.analysis.bam.bai'",
            "stageAs: 'control.analysis.bam'",
            "stageAs: 'control.analysis.bam.bai'",
        )),
    "broad calling uses target and matched IgG control":
        '-t "target.analysis.bam"' in broad
        and '-c "control.analysis.bam"' in broad,
    "broad calling preserves every NanoScope flag exactly":
        all(token in broad for token in (
            '-f BAMPE',
            '-n "sample"',
            '--llocal 100000',
            '--keep-dup 1',
            '--broad-cutoff 0.1',
            '--max-gap 1000',
            '--broad',
        )),
    "broad calling validates genome size before command interpolation":
        "validateMacsGenomeSize" in broad
        and "new BigInteger(text).signum() <= 0" in broad
        and "/[A-Za-z][A-Za-z0-9._-]*/" in broad
        and broad.index("validateMacsGenomeSize") < broad.index("macs2 callpeak"),
    "broad calling materializes valid empty output files after MACS2 succeeds":
        'touch "sample_peaks.broadPeak"' in broad
        and 'touch "sample_peaks.gappedPeak"' in broad
        and 'touch "sample_peaks.xls"' in broad,
    "broad outputs retain all raw MACS2 peak products":
        all(token in broad for token in (
            'sample_peaks.broadPeak',
            'sample_peaks.gappedPeak',
            'sample_peaks.xls',
        )),
    "broad outputs include metadata-keyed log and version records":
        broad.count("tuple val(meta)") >= 3
        and "macs2_broad.log" in broad
        and "macs2_broad_versions.yml" in broad,
    "MACS2_NARROW process is declared":
        "process MACS2_NARROW" in narrow,
    "narrow calling uses the pinned MACS2 environment and image":
        'conda "${projectDir}/envs/macs2.yml"' in narrow
        and f"container '{macs2_container}'" in narrow,
    "narrow calling remains matched-control BAMPE":
        '-t "target.analysis.bam"' in narrow
        and '-c "control.analysis.bam"' in narrow
        and '-f BAMPE' in narrow
        and '--llocal 100000' in narrow
        and '--keep-dup 1' in narrow,
    "narrow calling validates genome size and materializes empty outputs":
        "new BigInteger(text).signum() <= 0" in narrow
        and "/[A-Za-z][A-Za-z0-9._-]*/" in narrow
        and 'touch "sample_peaks.narrowPeak"' in narrow
        and 'touch "sample_summits.bed"' in narrow
        and 'touch "sample_peaks.xls"' in narrow,
    "narrow calling is explicitly optional":
        "val narrow_enabled" in narrow
        and "when:" in narrow
        and "narrow_enabled == true" in narrow,
    "narrow outputs are motif-labeled peaks, summits, log, and versions":
        all(token in narrow for token in (
            "motif_peaks",
            "motif_summits",
            "sample_peaks.narrowPeak",
            "sample_summits.bed",
            "macs2_narrow.log",
            "macs2_narrow_versions.yml",
        )),
    "FILTER_BLACKLIST process is declared":
        "process FILTER_BLACKLIST" in blacklist,
    "blacklist filtering uses the pinned BEDTools environment and image":
        'conda "${projectDir}/envs/bedtools.yml"' in blacklist
        and f"container '{bedtools_container}'" in blacklist,
    "blacklist input accepts the collected optional path without invalid arity":
        "arity: '0..1'" not in blacklist
        and "stageAs: 'blacklist/regions*.bed'" in blacklist,
    "blacklist filtering uses fixed paths and inverse intersection":
        'bedtools intersect' in blacklist
        and '-v' in blacklist
        and '-a "raw.broadPeak"' in blacklist
        and '-b "blacklist/regions.bed"' in blacklist,
    "absent blacklist copies raw peaks to a distinct final output":
        'cp "raw.broadPeak" "final.broadPeak"' in blacklist,
    "blacklist outputs include metadata-keyed final peaks, log, and versions":
        blacklist.count("tuple val(meta)") >= 3
        and "final.broadPeak" in blacklist
        and "filter_blacklist.log" in blacklist
        and "filter_blacklist_versions.yml" in blacklist,
    "peak workflow includes all three processes":
        all(name in peaks for name in (
            "MACS2_BROAD",
            "MACS2_NARROW",
            "FILTER_BLACKLIST",
        )),
    "peak workflow validates path-safe identifiers":
        "SAFE_ID" in peaks
        and "sample_id" in peaks
        and "control_id" in peaks,
    "only non-controls become target peak calls":
        "isControlMeta" in peaks
        and "!isControlMeta(meta)" in peaks,
    "controls are collected once into a reusable keyed map":
        "collect(flat: false)" in peaks
        and "controlsBySampleId" in peaks,
    "targets are paired to analysis controls by control_id":
        ".combine(controls_by_id)" in peaks
        and "meta.control_id" in peaks,
    "unmatched target controls fail with a clear target and control message":
        "no analysis BAM for control_id" in peaks
        and "referenced by target" in peaks,
    "self-control calls are explicitly rejected":
        "cannot use itself as an IgG control" in peaks,
    "an empty blacklist queue is collected into a reusable value":
        "blacklist.collect(flat: false)" in peaks,
    "broad peaks always flow through optional blacklist finalization":
        "FILTER_BLACKLIST(MACS2_BROAD.out.peaks" in peaks,
    "narrow calls receive the exact validated enable value":
        "MACS2_NARROW(paired_bams, safe_genome_size, narrow_enabled)" in peaks,
    "workflow emissions distinguish raw, final, and motif-only peaks":
        all(token in peaks for token in (
            "raw_broad_peaks =",
            "final_broad_peaks =",
            "motif_narrow_peaks =",
            "motif_summits =",
            "logs =",
            "versions =",
        )),
}

failed = [description for description, passed in checks.items() if not passed]
if failed:
    raise SystemExit(
        "FAIL: peak structural assertions failed:\n  - "
        + "\n  - ".join(failed)
    )

for source_name, source in (("macs2_broad.nf", broad), ("macs2_narrow.nf", narrow)):
    script = source.split("script:", 1)[1]
    if "${meta." in script or "${control_meta." in script:
        raise SystemExit(
            f"FAIL: {source_name} interpolates manifest metadata into its shell script"
        )

with (root / "assets/samples.example.csv").open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))
controls = {row["sample_id"]: row for row in rows if row["is_control"] == "true"}
targets = [row for row in rows if row["is_control"] == "false"]
expected_control = {
    "25K": "NX702_IgG",
    "50K": "NX703_IgG",
    "100K": "NX701_IgG",
}
assert len(targets) == 9
assert {row["input_group"] for row in targets} == set(expected_control)
for target in targets:
    control_id = target["control_id"]
    assert control_id == expected_control[target["input_group"]]
    assert control_id in controls
    assert controls[control_id]["assay_target"] == "IgG"
    assert controls[control_id]["input_group"] == target["input_group"]
PY

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/nanocut-peaks.XXXXXX")
trap 'rm -rf -- "${tmp_dir:?}"' EXIT
mkdir -p "$tmp_dir/fakebin" "$tmp_dir/direct"

cat > "$tmp_dir/fakebin/macs2" <<'PY'
#!/usr/bin/env python3
import os
import pathlib
import sys

args = sys.argv[1:]
if args == ["--version"]:
    print("macs2 2.2.9.1")
    raise SystemExit(0)
if not args or args[0] != "callpeak":
    raise SystemExit("fake macs2 only supports callpeak")
required = ("-t", "-c", "-g", "-f", "-n", "--llocal", "--keep-dup")
for option in required:
    if option not in args:
        raise SystemExit(f"missing required MACS2 option {option}")
if args[args.index("-t") + 1] == args[args.index("-c") + 1]:
    raise SystemExit("target and control BAM paths must differ")
name = args[args.index("-n") + 1]
log_path = os.environ.get("FAKE_MACS_ARGS")
if log_path:
    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write("\t".join(args) + "\n")
empty = os.environ.get("FAKE_EMPTY_PEAKS") == "1"
if "--broad" in args:
    peak_text = "" if empty else "chrMini\t10\t30\tpeak\t50\t.\t1\t2\t3\n"
    pathlib.Path(f"{name}_peaks.broadPeak").write_text(peak_text)
    pathlib.Path(f"{name}_peaks.gappedPeak").write_text(peak_text)
else:
    peak_text = "" if empty else "chrMini\t12\t28\tpeak\t60\t.\t2\t3\t4\t8\n"
    pathlib.Path(f"{name}_peaks.narrowPeak").write_text(peak_text)
    pathlib.Path(f"{name}_summits.bed").write_text(
        "" if empty else "chrMini\t20\t21\tpeak\t4\n"
    )
pathlib.Path(f"{name}_peaks.xls").write_text("# fake MACS2 output\n")
PY
chmod +x "$tmp_dir/fakebin/macs2"

cat > "$tmp_dir/fakebin/bedtools" <<'PY'
#!/usr/bin/env python3
import pathlib
import sys

args = sys.argv[1:]
if args == ["--version"]:
    print("bedtools v2.31.1")
    raise SystemExit(0)
if not args or args[0] != "intersect" or "-v" not in args:
    raise SystemExit("fake bedtools only supports intersect -v")
a_path = pathlib.Path(args[args.index("-a") + 1])
b_path = pathlib.Path(args[args.index("-b") + 1])
excluded = []
for line in b_path.read_text().splitlines():
    if line.strip():
        chrom, start, end, *_ = line.split("\t")
        excluded.append((chrom, int(start), int(end)))
for line in a_path.read_text().splitlines():
    if not line.strip():
        continue
    chrom, start, end, *_ = line.split("\t")
    start_i, end_i = int(start), int(end)
    overlaps = any(
        chrom == b_chrom and start_i < b_end and b_start < end_i
        for b_chrom, b_start, b_end in excluded
    )
    if not overlaps:
        print(line)
PY
chmod +x "$tmp_dir/fakebin/bedtools"

: > "$tmp_dir/direct/target.bam"
: > "$tmp_dir/direct/control.bam"
printf 'chrMini\t0\t100\n' > "$tmp_dir/direct/blacklist.bed"

(
  cd "$tmp_dir/direct"
  PATH="$tmp_dir/fakebin:$PATH" \
  FAKE_MACS_ARGS="$tmp_dir/macs.args" \
  FAKE_EMPTY_PEAKS=1 \
    macs2 callpeak \
      -t "target.bam" \
      -c "control.bam" \
      -g "hs" \
      -f BAMPE \
      -n "sample" \
      --llocal 100000 \
      --keep-dup 1 \
      --broad-cutoff 0.1 \
      --max-gap 1000 \
      --broad
  PATH="$tmp_dir/fakebin:$PATH" \
    bedtools intersect \
      -v \
      -a "sample_peaks.broadPeak" \
      -b "blacklist.bed" \
      > "final.broadPeak"
)

[[ -f "$tmp_dir/direct/sample_peaks.broadPeak" ]]
[[ ! -s "$tmp_dir/direct/sample_peaks.broadPeak" ]]
[[ -f "$tmp_dir/direct/final.broadPeak" ]]
[[ ! -s "$tmp_dir/direct/final.broadPeak" ]]
grep -F $'callpeak\t-t\ttarget.bam\t-c\tcontrol.bam\t-g\ths\t-f\tBAMPE' \
  "$tmp_dir/macs.args" >/dev/null

if ! command -v nextflow >/dev/null 2>&1; then
  printf '%s\n' \
    'SKIP: Nextflow runtime unavailable; strong static checks and direct fake-tool empty-peak checks passed, but no DSL2 execution was run.'
  exit 0
fi

mkdir -p "$tmp_dir/runtime/input"
for sample in IgG TARGET; do
  : > "$tmp_dir/runtime/input/${sample}.bam"
  : > "$tmp_dir/runtime/input/${sample}.bam.bai"
done

cat > "$tmp_dir/runtime/main.nf" <<EOF
nextflow.enable.dsl = 2

include { PEAKS } from '${repo_root}/subworkflows/local/peaks'

workflow {
    analysis_bams = Channel.of(
        tuple(
            [sample_id: 'IgG', is_control: true, control_id: null],
            file('${tmp_dir}/runtime/input/IgG.bam'),
            file('${tmp_dir}/runtime/input/IgG.bam.bai')
        ),
        tuple(
            [sample_id: 'TARGET', is_control: false, control_id: 'IgG'],
            file('${tmp_dir}/runtime/input/TARGET.bam'),
            file('${tmp_dir}/runtime/input/TARGET.bam.bai')
        )
    )
    PEAKS(analysis_bams, Channel.empty(), 'hs', true)
}
EOF

cat > "$tmp_dir/runtime/nextflow.config" <<EOF
params.outdir = '${tmp_dir}/runtime/results'
process.executor = 'local'
process.cpus = 1
process.memory = '1 GB'
process.time = '10m'
conda.enabled = false
docker.enabled = false
EOF

(
  cd "$tmp_dir/runtime"
  PATH="$tmp_dir/fakebin:$PATH" \
  FAKE_MACS_ARGS="$tmp_dir/runtime.macs.args" \
  FAKE_EMPTY_PEAKS=1 \
    nextflow run main.nf -c nextflow.config -ansi-log false \
      -with-trace "$tmp_dir/runtime/trace.txt"
)

[[ -f "$tmp_dir/runtime/results/peaks/TARGET/broad/raw/sample_peaks.broadPeak" ]]
[[ -f "$tmp_dir/runtime/results/peaks/TARGET/broad/final/final.broadPeak" ]]
[[ -f "$tmp_dir/runtime/results/peaks/TARGET/narrow_motif_qc/sample_peaks.narrowPeak" ]]
[[ ! -s "$tmp_dir/runtime/results/peaks/TARGET/broad/raw/sample_peaks.broadPeak" ]]
[[ ! -s "$tmp_dir/runtime/results/peaks/TARGET/broad/final/final.broadPeak" ]]
[[ $(grep -c 'MACS2_BROAD' "$tmp_dir/runtime/trace.txt") -eq 1 ]]
[[ $(grep -c 'MACS2_NARROW' "$tmp_dir/runtime/trace.txt") -eq 1 ]]
[[ $(grep -c 'FILTER_BLACKLIST' "$tmp_dir/runtime/trace.txt") -eq 1 ]]
[[ $(wc -l < "$tmp_dir/runtime.macs.args") -eq 2 ]]

cat > "$tmp_dir/runtime/narrow_disabled.nf" <<EOF
nextflow.enable.dsl = 2

include { PEAKS } from '${repo_root}/subworkflows/local/peaks'

workflow {
    analysis_bams = Channel.of(
        tuple(
            [sample_id: 'IgG', assay_target: 'IgG', is_control: true, control_id: null],
            file('${tmp_dir}/runtime/input/IgG.bam'),
            file('${tmp_dir}/runtime/input/IgG.bam.bai')
        ),
        tuple(
            [sample_id: 'TARGET', assay_target: 'CTCF', is_control: false, control_id: 'IgG'],
            file('${tmp_dir}/runtime/input/TARGET.bam'),
            file('${tmp_dir}/runtime/input/TARGET.bam.bai')
        )
    )
    PEAKS(analysis_bams, Channel.empty(), '1000', false)
}
EOF

rm -rf "$tmp_dir/runtime/results"
(
  cd "$tmp_dir/runtime"
  PATH="$tmp_dir/fakebin:$PATH" \
  FAKE_MACS_ARGS="$tmp_dir/runtime-disabled.macs.args" \
    nextflow run narrow_disabled.nf -c nextflow.config -ansi-log false \
      -with-trace "$tmp_dir/runtime/disabled-trace.txt"
)

[[ -f "$tmp_dir/runtime/results/peaks/TARGET/broad/raw/sample_peaks.broadPeak" ]]
[[ -f "$tmp_dir/runtime/results/peaks/TARGET/broad/final/final.broadPeak" ]]
[[ ! -e "$tmp_dir/runtime/results/peaks/TARGET/narrow_motif_qc/sample_peaks.narrowPeak" ]]
[[ $(wc -l < "$tmp_dir/runtime-disabled.macs.args") -eq 1 ]]

cat > "$tmp_dir/runtime/unmatched.nf" <<EOF
nextflow.enable.dsl = 2

include { PEAKS } from '${repo_root}/subworkflows/local/peaks'

workflow {
    analysis_bams = Channel.of(
        tuple(
            [sample_id: 'TARGET', is_control: false, control_id: 'MISSING'],
            file('${tmp_dir}/runtime/input/TARGET.bam'),
            file('${tmp_dir}/runtime/input/TARGET.bam.bai')
        )
    )
    PEAKS(analysis_bams, Channel.empty(), '1000', false)
}
EOF

if (
  cd "$tmp_dir/runtime"
  PATH="$tmp_dir/fakebin:$PATH" \
    nextflow run unmatched.nf -c nextflow.config -ansi-log false \
      > "$tmp_dir/unmatched.log" 2>&1
); then
  printf 'FAIL: unmatched control_id unexpectedly succeeded\n' >&2
  exit 1
fi
grep -F \
  'no analysis BAM for control_id MISSING referenced by target TARGET' \
  "$tmp_dir/unmatched.log" >/dev/null

printf '%s\n' \
  'PASS: matched-control broad/optional narrow peak calls, empty outputs, and absent blacklist handling executed successfully.'
