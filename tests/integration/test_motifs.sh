#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"

required_files=(
  modules/local/prepare_motif_sequences.nf
  modules/local/ame.nf
  modules/local/streme.nf
  modules/local/fimo.nf
  modules/local/motif_summary.nf
  subworkflows/local/motifs.nf
)

for required_file in "${required_files[@]}"; do
  if [[ ! -f "$required_file" ]]; then
    printf 'FAIL: missing %s\n' "$required_file" >&2
    exit 1
  fi
done

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/nanocut-motifs.XXXXXX")
trap 'rm -rf -- "${tmp_dir:?}"' EXIT

python3 - "$tmp_dir" <<'PY'
import json
import re
import subprocess
import sys
from pathlib import Path

root = Path.cwd()
tmp = Path(sys.argv[1])
prepare = (root / "modules/local/prepare_motif_sequences.nf").read_text()
ame = (root / "modules/local/ame.nf").read_text()
streme = (root / "modules/local/streme.nf").read_text()
fimo = (root / "modules/local/fimo.nf").read_text()
summary = (root / "modules/local/motif_summary.nf").read_text()
motifs = (root / "subworkflows/local/motifs.nf").read_text()

meme_container = "quay.io/biocontainers/meme:5.5.7--pl5321h1ca524f_3"
bedtools_container = "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"
script_start = prepare.index('    """\n    set -euo pipefail')
script_end = prepare.rindex('    """')
interpolated_script = prepare[script_start:script_end]

checks = {
    "PREPARE_MOTIF_SEQUENCES process is declared":
        "process PREPARE_MOTIF_SEQUENCES" in prepare,
    "sequence preparation uses the pinned BEDTools environment and image":
        'conda "${projectDir}/envs/bedtools.yml"' in prepare
        and f"container '{bedtools_container}'" in prepare,
    "optional blacklist staging does not declare zero-based path arity":
        "arity: '0..1'" not in prepare
        and "stageAs: 'blacklist/regions*.bed'" in prepare,
    "Groovy-interpolated shell escapes every literal dollar sign":
        re.search(r"(?<!\\)\$(?!\{)", interpolated_script) is None,
    "sequence preparation validates a positive integral total window":
        "validateMotifWindow" in prepare
        and "signum() <= 0" in prepare
        and "motif_window must be a positive integer" in prepare,
    "summit and broad-midpoint anchor modes are explicit":
        "summit" in prepare and "midpoint" in prepare
        and "unsupported motif anchor mode" in prepare,
    "narrow summits are restricted to final non-blacklisted peak regions":
        'bedtools intersect' in prepare
        and '-u -a "motif.anchor.bed" -b "final.peaks.bed"' in prepare
        and '-v -a "anchors.in.final.peaks.bed" -b "blacklist/regions.bed"' in prepare
        and '"usable.anchors.bed"' in prepare,
    "reference bounds are calculated from the staged FASTA":
        'reference.genome' in prepare
        and 'reference.fa' in prepare
        and 'make_windows.awk' in prepare,
    "foreground windows and FASTA use fixed staged paths":
        'foreground.windows.bed' in prepare
        and 'foreground.fa' in prepare
        and 'bedtools getfasta' in prepare,
    "background excludes peaks and the optional blacklist":
        'background_exclusions.bed' in prepare
        and '-excl "background_exclusions.bed"' in prepare
        and 'blacklist/regions.bed' in prepare,
    "background shuffling is reproducible and mutually non-overlapping":
        "-seed 1729" in prepare
        and "-noOverlapping" in prepare
        and "-maxTries 10000" in prepare,
    "background candidates preserve lengths and are approximately GC matched":
        'candidate.features.tsv' in prepare
        and 'foreground.features.tsv' in prepare
        and 'match_gc.awk' in prepare
        and "candidate_length[candidate_index] != foreground_length[foreground_index]" in prepare,
    "empty foregrounds produce a recorded no_peaks result":
        "no_peaks" in prepare
        and 'motif_sequence_status.tsv' in prepare,
    "AME process is declared":
        "process AME" in ame,
    "AME uses the pinned MEME environment and image":
        'conda "${projectDir}/envs/meme.yml"' in ame
        and f"container '{meme_container}'" in ame,
    "AME uses a separate shuffled control and requests the complete database":
        'ame --oc "ame"' in ame
        and '--control "background.fa"' in ame
        and "--method fisher" in ame
        and "--scoring avg" in ame
        and "--evalue-report-threshold 1e300" in ame
        and "--pvalue-report-threshold" not in ame,
    "AME complete-database provenance is added only after a successful run":
        "# motif_qc_complete_database=true" in ame
        and ame.index('ame --oc "ame"') < ame.index("# motif_qc_complete_database=true"),
    "AME records empty peaks without invoking the tool":
        'if [[ ! -s "foreground.fa" ]]' in ame
        and 'status\\\\tno_peaks' in ame,
    "STREME process is declared":
        "process STREME" in streme,
    "STREME uses the pinned MEME environment and modern positive/negative CLI":
        'conda "${projectDir}/envs/meme.yml"' in streme
        and f"container '{meme_container}'" in streme
        and 'streme --oc "streme"' in streme
        and '--p "foreground.fa"' in streme
        and '--n "background.fa"' in streme,
    "FIMO process is declared":
        "process FIMO" in fimo,
    "FIMO uses the pinned MEME environment and modern CLI":
        'conda "${projectDir}/envs/meme.yml"' in fimo
        and f"container '{meme_container}'" in fimo
        and 'fimo --oc "fimo"' in fimo
        and '"motifs.meme" "foreground.fa"' in fimo,
    "all MEME wrappers retain result directories, status, logs, and versions":
        all(
            all(token in source for token in ("emit: results", "emit: status", "emit: logs", "emit: versions"))
            for source in (ame, streme, fimo)
        ),
    "MEME result directories publish exactly at motifs/sample/tool":
        all(
            'publishDir "${params.outdir}/motifs/${meta.sample_id}"' in source
            and "saveAs:" not in source
            for source in (ame, streme, fimo)
        ),
    "MOTIF_SUMMARY process is declared":
        "process MOTIF_SUMMARY" in summary,
    "motif summary wraps the tested repository utility":
        'conda "${projectDir}/envs/python.yml"' in summary
        and "container 'python:3.12.3-slim-bookworm'" in summary
        and "motif_qc.py" in summary
        and '--ame "ame/ame.tsv"' in summary
        and '--fimo "fimo/fimo.tsv"' in summary,
    "expected-motif regex is transported without shell interpolation":
        "java.util.Base64.getEncoder().encodeToString" in summary
        and 'expected_motif.b64' in summary
        and '${meta.expected_motif}' not in summary,
    "summary emits duplicate-safe metadata-keyed JSON TSV and positions":
        summary.count("tuple val(meta)") >= 2
        and "expected_motif_qc.json" in summary
        and "expected_motif_qc.tsv" in summary
        and "expected_motif_qc.motif_hit_positions.tsv" in summary,
    "MOTIFS includes all five modules":
        all(name in motifs for name in (
            "PREPARE_MOTIF_SEQUENCES", "AME", "STREME", "FIMO", "MOTIF_SUMMARY",
        )),
    "motif workflow validates safe target metadata and expected-motif regex":
        "validateMotifMeta" in motifs
        and "is_control must be false" in motifs
        and "expected_motif must not be blank" in motifs
        and "Pattern.compile" in motifs,
    "motif workflow rejects duplicate broad peaks and summits":
        "duplicate ${label} sample_id" in motifs
        and "uniqueMotifRows(rows, 'final motif peak')" in motifs
        and "uniqueMotifRows(rows, 'motif summit')" in motifs,
    "narrow summits are selected when enabled and broad midpoints otherwise":
        "motif_use_narrow_peaks" in motifs
        and "'summit'" in motifs
        and "'midpoint'" in motifs
        and "missing motif summit" in motifs,
    "collected peak inputs combine as maps rather than flattenable row lists":
        "safe_broad_by_sample" in motifs
        and "safe_summits_by_sample" in motifs
        and "broadBySample.values().collect" in motifs,
    "AME and FIMO results are joined by sample with strict duplicate/mismatch handling":
        "failOnDuplicate: true" in motifs
        and "failOnMismatch: true" in motifs
        and "motif result metadata mismatch" in motifs,
    "Task 9 receives exact metadata-keyed motif metric tuples":
        "motif_metrics =" in motifs
        and "tuple(meta, metricsTsv)" in motifs,
    "workflow emits all retained motif products and versions":
        all(token in motifs for token in (
            "known_motifs =", "de_novo_motifs =", "scans =",
            "expected_motif_qc =", "motif_metrics =", "versions =",
        )),
}

failed = [description for description, passed in checks.items() if not passed]
if failed:
    raise SystemExit(
        "FAIL: motif structural assertions failed:\n  - "
        + "\n  - ".join(failed)
    )

for source_name, source in (
    ("prepare_motif_sequences.nf", prepare),
    ("ame.nf", ame),
    ("streme.nf", streme),
    ("fimo.nf", fimo),
    ("motif_summary.nf", summary),
):
    shell = source.split("script:", 1)[1]
    if "${meta." in shell:
        raise SystemExit(
            f"FAIL: {source_name} interpolates manifest metadata into its shell script"
        )

window_match = re.search(
    r'cat > "make_windows\.awk" <<\'AWK\'\n(.*?)\nAWK',
    prepare,
    re.DOTALL,
)
gc_match = re.search(
    r'cat > "match_gc\.awk" <<\'AWK\'\n(.*?)\nAWK',
    prepare,
    re.DOTALL,
)
bed3_match = re.search(
    r'cat > "normalize_bed3\.awk" <<\'AWK\'\n(.*?)\nAWK',
    prepare,
    re.DOTALL,
)
if not window_match or not gc_match or not bed3_match:
    raise SystemExit(
        "FAIL: could not extract motif window/GC/BED3 normalization programs"
    )

def render_nextflow_gstring(text):
    return (
        text
        .replace("\\$", "$")
        .replace("\\\\t", "\\t")
        .replace("\\\\n", "\\n")
    )

(tmp / "make_windows.awk").write_text(
    render_nextflow_gstring(window_match.group(1)) + "\n",
    encoding="utf-8",
)
(tmp / "match_gc.awk").write_text(
    render_nextflow_gstring(gc_match.group(1)) + "\n",
    encoding="utf-8",
)
(tmp / "normalize_bed3.awk").write_text(
    render_nextflow_gstring(bed3_match.group(1)) + "\n",
    encoding="utf-8",
)

(tmp / "mixed_width_exclusions.bed").write_text(
    "chr1\t10115\t10430\tsample_peak_57\t80\t.\t2.9552\t11.4229\t8.0468\n"
    "# blacklist regions use BED3\n"
    "chr1\t20000\t20200\n",
    encoding="utf-8",
)
bed3_result = subprocess.run(
    [
        "awk", "-f", str(tmp / "normalize_bed3.awk"),
        str(tmp / "mixed_width_exclusions.bed"),
    ],
    capture_output=True,
    text=True,
    check=True,
)
assert bed3_result.stdout.splitlines() == [
    "chr1\t10115\t10430",
    "chr1\t20000\t20200",
]

(tmp / "reference.genome").write_text("chrMini\t500\nchrShort\t80\n", encoding="utf-8")
(tmp / "summits.bed").write_text(
    "chrMini\t0\t1\n"
    "chrMini\t250\t251\n"
    "chrMini\t499\t500\n"
    "chrShort\t30\t31\n",
    encoding="utf-8",
)
window_result = subprocess.run(
    [
        "awk", "-v", "mode=summit", "-v", "window=200",
        "-f", str(tmp / "make_windows.awk"),
        str(tmp / "reference.genome"), str(tmp / "summits.bed"),
    ],
    capture_output=True,
    text=True,
    check=True,
)
assert window_result.stdout.splitlines() == [
    "chrMini\t0\t200\tpeak_000001",
    "chrMini\t150\t350\tpeak_000002",
    "chrMini\t300\t500\tpeak_000003",
    "chrShort\t0\t80\tpeak_000004",
]
for line in window_result.stdout.splitlines()[:3]:
    chrom, start, end, _ = line.split("\t")
    assert chrom == "chrMini"
    assert int(end) - int(start) == 200

(tmp / "broad.bed").write_text("chrMini\t0\t20\n", encoding="utf-8")
midpoint_result = subprocess.run(
    [
        "awk", "-v", "mode=midpoint", "-v", "window=200",
        "-f", str(tmp / "make_windows.awk"),
        str(tmp / "reference.genome"), str(tmp / "broad.bed"),
    ],
    capture_output=True,
    text=True,
    check=True,
)
assert midpoint_result.stdout.strip() == "chrMini\t0\t200\tpeak_000001"

(tmp / "foreground.features.tsv").write_text(
    "peak_000001\t200\t0.20\n"
    "peak_000002\t200\t0.75\n",
    encoding="utf-8",
)
(tmp / "candidate.features.tsv").write_text(
    "chrMini\t0\t200\tcandidate_1\t200\t0.74\n"
    "chrMini\t200\t400\tcandidate_2\t200\t0.18\n"
    "chrMini\t300\t380\tcandidate_wrong_length\t80\t0.75\n",
    encoding="utf-8",
)
gc_result = subprocess.run(
    [
        "awk", "-f", str(tmp / "match_gc.awk"),
        str(tmp / "foreground.features.tsv"),
        str(tmp / "candidate.features.tsv"),
    ],
    capture_output=True,
    text=True,
    check=True,
)
assert gc_result.stdout.splitlines() == [
    "chrMini\t200\t400\tpeak_000001",
    "chrMini\t0\t200\tpeak_000002",
]

complete_ame = tmp / "complete.ame.tsv"
complete_ame.write_text(
    "# motif_qc_complete_database=true\n"
    "rank\tmotif_ID\tmotif_Alt_ID\tE-value\tpos\n"
    "1\tMA0001.1\tOTHER\t0.9\t2\n",
    encoding="utf-8",
)
empty_fimo = tmp / "empty.fimo.tsv"
empty_fimo.write_text(
    "# motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tp-value\tq-value\n",
    encoding="utf-8",
)
complete_prefix = tmp / "complete"
subprocess.run(
    [
        sys.executable, str(root / "bin" / "motif_qc.py"),
        "--ame", str(complete_ame),
        "--fimo", str(empty_fimo),
        "--expected-motif", "^CTCF$",
        "--window", "200",
        "--output-prefix", str(complete_prefix),
    ],
    check=True,
)
complete_payload = json.loads((tmp / "complete.json").read_text(encoding="utf-8"))
assert complete_payload["status"] == "motif_not_found"
assert complete_payload["ame_database_complete"] is True

empty_ame = tmp / "empty.ame.tsv"
empty_ame.write_text(
    "rank\tmotif_ID\tmotif_Alt_ID\tp-value\tE-value\tpos\tneg\n"
    "1\t__NO_PEAKS__\t__NO_PEAKS__\t1\t1\t0\t0\n",
    encoding="utf-8",
)
empty_prefix = tmp / "empty"
subprocess.run(
    [
        sys.executable, str(root / "bin" / "motif_qc.py"),
        "--ame", str(empty_ame),
        "--fimo", str(empty_fimo),
        "--expected-motif", "^CTCF$",
        "--window", "200",
        "--output-prefix", str(empty_prefix),
    ],
    check=True,
)
empty_payload = json.loads((tmp / "empty.json").read_text(encoding="utf-8"))
assert empty_payload["status"] == "no_peaks"
PY

if ! command -v nextflow >/dev/null 2>&1; then
  printf '%s\n' \
    'SKIP: Nextflow runtime unavailable; strong static checks and direct bounded-window/GC-match/provenance/empty-peak fixtures passed, but no DSL2 execution was run.'
  exit 0
fi

mkdir -p "$tmp_dir/runtime"
cat > "$tmp_dir/runtime/reference.fa" <<'EOF'
>chrMini
ACGTACGTACGTACGTACGT
EOF
cat > "$tmp_dir/runtime/motifs.meme" <<'EOF'
MEME version 5

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF MA0001.1 TEST
letter-probability matrix: alength= 4 w= 4 nsites= 2 E= 0
0.9 0.03 0.03 0.04
0.03 0.9 0.03 0.04
0.03 0.03 0.9 0.04
0.04 0.03 0.03 0.9
EOF
cat > "$tmp_dir/runtime/main.nf" <<EOF
nextflow.enable.dsl = 2

include { MOTIFS } from '${repo_root}/subworkflows/local/motifs'

workflow {
    MOTIFS(
        Channel.empty(),
        Channel.empty(),
        Channel.of(file('${tmp_dir}/runtime/reference.fa')),
        Channel.empty(),
        Channel.of(file('${tmp_dir}/runtime/motifs.meme')),
        Channel.value(200),
        Channel.value(false)
    )
}
EOF
cat > "$tmp_dir/runtime/nextflow.config" <<EOF
params.outdir = '${tmp_dir}/runtime/results'
process.executor = 'local'
conda.enabled = false
docker.enabled = false
EOF
(
  cd "$tmp_dir/runtime"
  nextflow run main.nf -c nextflow.config -ansi-log false
)

printf '%s\n' \
  'PASS: motif DSL2 parsed/executed with empty target channels and direct bounded-window/GC-match/provenance/empty-peak fixtures passed.'
