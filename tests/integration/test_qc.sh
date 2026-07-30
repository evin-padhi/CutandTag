#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"

required_files=(
  modules/local/bam_to_fragments.nf
  modules/local/filtered_bam_qc.nf
  modules/local/peak_qc.nf
  modules/local/tss_enrichment.nf
  modules/local/multiqc.nf
  subworkflows/local/qc.nf
)

for required_file in "${required_files[@]}"; do
  if [[ ! -f "$required_file" ]]; then
    printf 'FAIL: missing %s\n' "$required_file" >&2
    exit 1
  fi
done

python3 - <<'PY'
from pathlib import Path
import re

root = Path.cwd()
bam_to_fragments = (root / "modules/local/bam_to_fragments.nf").read_text()
filtered_bam_qc = (root / "modules/local/filtered_bam_qc.nf").read_text()
peak_qc = (root / "modules/local/peak_qc.nf").read_text()
tss = (root / "modules/local/tss_enrichment.nf").read_text()
multiqc = (root / "modules/local/multiqc.nf").read_text()
qc = (root / "subworkflows/local/qc.nf").read_text()

samtools_image = "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
deeptools_image = "quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0"
multiqc_image = "quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0"

checks = {
    "BAM_TO_FRAGMENTS process is declared":
        "process BAM_TO_FRAGMENTS" in bam_to_fragments,
    "fragment conversion uses the pinned SAMtools environment and image":
        'conda "${projectDir}/envs/samtools.yml"' in bam_to_fragments
        and f"container '{samtools_image}'" in bam_to_fragments,
    "fragment conversion consumes only the filtered BAM contract":
        "filtered_bam" in bam_to_fragments
        and "stageAs: 'qc.filtered.bam'" in bam_to_fragments
        and "analysis_bam" not in bam_to_fragments,
    "fragment conversion reasserts proper-pair primary alignment filters":
        "samtools view" in bam_to_fragments
        and "-f 2" in bam_to_fragments
        and "-F 2304" in bam_to_fragments,
    "fragment conversion is name-collated then deterministically sorted":
        "samtools sort" in bam_to_fragments
        and "-n" in bam_to_fragments
        and "LC_ALL=C sort" in bam_to_fragments,
    "fragment conversion emits one required name-bearing BEDPE record per pair":
        "fragment_pairs.awk" in bam_to_fragments
        and "expected exactly two primary proper-pair alignments" in bam_to_fragments
        and 'print c1, s1, e1, c2, s2, e2, current' in bam_to_fragments,
    "fragment conversion validates exactly one first and one second mate":
        "has_flag" in bam_to_fragments
        and "has_flag(\\$2, 64)" in bam_to_fragments
        and "has_flag(\\$2, 128)" in bam_to_fragments
        and "first_bit == second_bit" in bam_to_fragments
        and "each alignment must set exactly one of 0x40 and 0x80" in bam_to_fragments
        and "expected one 0x40 and one 0x80 mate" in bam_to_fragments,
    "fragment conversion emits metadata-keyed fragments and versions":
        bam_to_fragments.count("tuple val(meta)") >= 2
        and "fragments.bedpe" in bam_to_fragments
        and "versions" in bam_to_fragments,
    "FILTERED_BAM_QC process is declared":
        "process FILTERED_BAM_QC" in filtered_bam_qc,
    "filtered BAM QC uses the pinned SAMtools environment and image":
        'conda "${projectDir}/envs/samtools.yml"' in filtered_bam_qc
        and f"container '{samtools_image}'" in filtered_bam_qc,
    "filtered BAM QC counts first-only and second-only identities":
        'samtools view -c "qc.filtered.bam"' in filtered_bam_qc
        and 'samtools view -c -f 64 -F 128 "qc.filtered.bam"' in filtered_bam_qc
        and 'samtools view -c -f 128 -F 64 "qc.filtered.bam"' in filtered_bam_qc
        and "filtered_flags.awk" in filtered_bam_qc,
    "filtered BAM QC validates complete pairs and emits keyed counts":
        "filtered first/second mate counts differ" in filtered_bam_qc
        and "first-only plus second-only counts do not equal total records" in filtered_bam_qc
        and "filtered_bam_qc.tsv" in filtered_bam_qc
        and filtered_bam_qc.count("tuple val(meta)") >= 2,
    "PEAK_QC process is declared":
        "process PEAK_QC" in peak_qc,
    "peak QC uses the project Python environment and pinned Python image":
        'conda "${projectDir}/envs/python.yml"' in peak_qc
        and "container 'python:3.12.3-slim-bookworm'" in peak_qc,
    "peak QC consumes name-bearing fragments and final broad peaks":
        "stageAs: 'fragments.bedpe'" in peak_qc
        and "stageAs: 'final.broadPeak'" in peak_qc,
    "peak QC wraps the tested repository CLI":
        "peak_qc.py" in peak_qc
        and '--fragments "fragments.bedpe"' in peak_qc
        and '--peaks "final.broadPeak"' in peak_qc,
    "peak QC retains JSON, metric, width, and fragments-per-peak outputs":
        all(token in peak_qc for token in (
            ".peak_qc.json",
            ".peak_qc.tsv",
            ".peak_qc.width_histogram.tsv",
            ".peak_qc.fragments_per_peak.tsv",
        )),
    "TSS_ENRICHMENT process is declared":
        "process TSS_ENRICHMENT" in tss,
    "TSS enrichment uses the pinned deepTools environment and image":
        'conda "${projectDir}/envs/deeptools.yml"' in tss
        and f"container '{deeptools_image}'" in tss,
    "GTF conversion selects strand-aware transcript TSS coordinates":
        '$3 == "transcript"' in tss
        and '$7 == "+"' in tss
        and '$7 == "-"' in tss
        and "start = $4 - 1" in tss
        and "start = $5 - 1" in tss,
    "precomputed TSS BED is validated as strand-aware BED6":
        "BED6" in tss
        and '$6 != "+" && $6 != "-"' in tss,
    "deepTools computes a reference-point TSS matrix and profile":
        "computeMatrix reference-point" in tss
        and "--referencePoint TSS" in tss
        and "plotProfile" in tss
        and '--outFileNameData "${outputStem}.tss_profile.tsv"' in tss
        and 'path("*.tss_profile.tsv")' in tss,
    "MULTIQC and custom-content formatter processes are declared":
        all(name in multiqc for name in (
            "process DEMUX_QC_CUSTOM",
            "process LIBRARY_QC_CUSTOM",
            "process MOTIF_QC_CUSTOM",
            "process MULTIQC",
        )),
    "library formatter consumes filtered counts and real markdup JSON keys":
        "filtered_bam_qc.tsv" in multiqc
        and '"DUPLICATE TOTAL"' in multiqc
        and '"EXAMINED"' in multiqc
        and '"READ"' in multiqc
        and '"ESTIMATED LIBRARY SIZE"' in multiqc
        and 'count_for("duplicates")' not in multiqc,
    "library formatter validates markdup numeric fields":
        "require_non_negative_number" in multiqc
        and "missing required {source} field" in multiqc
        and "must be a JSON number" in multiqc,
    "library formatter uses alignment stats and insert-size rows":
        '"raw total sequences"' in multiqc
        and "mapq_filtered_fraction" in multiqc
        and "insert_size_mean" in multiqc
        and "insert_size_median" in multiqc
        and "insert_size_distribution.tsv" in multiqc,
    "motif formatter carries metadata identity into long-form motif TSV":
        "motifQcSampleId" in multiqc
        and "stageAs: 'expected_motif_qc.tsv'" in multiqc
        and '"sample_id": sample_id' in multiqc
        and "duplicate motif sample_id" in multiqc,
    "MultiQC uses the pinned project environment and image":
        'conda "${projectDir}/envs/multiqc.yml"' in multiqc
        and f"container '{multiqc_image}'" in multiqc,
    "optional MultiQC collections do not declare zero-based path arity":
        "arity: '0..*'" not in multiqc
        and all(stage in multiqc for stage in (
            "stageAs: 'fastqc??/*'",
            "stageAs: 'demux??/*'",
            "stageAs: 'library??/*'",
            "stageAs: 'insert??/*'",
            "stageAs: 'peak_qc??/*'",
            "stageAs: 'motif??/*'",
            "stageAs: 'tss??/*'",
        )),
    "custom report content covers demux, FRiP, peak count, motif, and library QC":
        all(token in multiqc for token in (
            "nanocut_demultiplex",
            "assigned_fraction",
            "nanocut_peak_qc",
            "frip",
            "peak_count",
            "nanocut_motif_qc",
            "expected_motif_status",
            "nanocut_library_qc",
            "properly_paired_percent",
            "mapq_filtered_fraction",
            "nanocut_insert_size",
        )),
    "MultiQC writes a collision-safe combined target summary":
        "combined_target_qc.tsv" in multiqc
        and "peak_qc??/*" in multiqc
        and "duplicate peak-QC sample_id" in multiqc,
    "MultiQC retains HTML, data, custom content, and versions":
        "multiqc_report.html" in multiqc
        and "multiqc_data" in multiqc
        and "multiqc_custom_content" in multiqc
        and "multiqc_versions.yml" in multiqc,
    "QC subworkflow includes every Task 9 process":
        all(name in qc for name in (
            "BAM_TO_FRAGMENTS",
            "FILTERED_BAM_QC",
            "PEAK_QC",
            "TSS_ENRICHMENT",
            "DEMUX_QC_CUSTOM",
            "LIBRARY_QC_CUSTOM",
            "MOTIF_QC_CUSTOM",
            "MULTIQC",
        ))
        and (
            "{ DEMUX_QC_CUSTOM; LIBRARY_QC_CUSTOM; "
            "MOTIF_QC_CUSTOM; MULTIQC }"
        ) in qc,
    "QC validates control metadata and excludes IgG from default FRiP":
        "is_control must be a boolean" in qc
        and "!meta.is_control" in qc
        and "target_filtered_bams" in qc,
    "FRiP pairs filtered target BAMs with only final broad peaks":
        "safe_filtered_bams" in qc
        and "safe_final_broad_peaks" in qc
        and "failOnMismatch: true" in qc,
    "optional annotations are collected to reusable values without deadlock":
        qc.count("collect(flat: false)") >= 2
        and "skipped_no_annotation" in qc
        and "annotation_mode != 'none'" in qc,
    "TSS BED takes precedence over GTF when both are present":
        "tss_rows ? 'bed'" in qc,
    "TSS profile tuple retains the profile-data table":
        "meta, bed, matrix, matrixTable, profile, profileTable, status" in qc,
    "IgG and target library metrics both feed the report":
        "FILTERED_BAM_QC(safe_filtered_bams)" in qc
        and "LIBRARY_QC_CUSTOM(library_qc_inputs)" in qc,
    "library metrics join their own filtered-BAM counts by sample ID":
        "keyed_library_metrics" in qc
        and "keyed_filtered_bam_qc" in qc
        and "library QC metadata mismatch" in qc
        and "failOnMismatch: true" in qc,
    "Task 10 motif contract is a metadata-keyed TSV tuple":
        "MOTIF_QC_CUSTOM(safe_motif_metrics)" in qc
        and "meta, motifTsv" in qc
        and "Task 10 contract: tuple(meta, expected_motif_qc_tsv)" in qc,
    "QC emits target QC, combined summary, optional TSS, report, and versions":
        all(token in qc for token in (
            "target_qc =",
            "combined_summary =",
            "tss_profiles =",
            "multiqc_report =",
            "multiqc_data =",
            "versions =",
        )),
}

failed = [description for description, passed in checks.items() if not passed]
if failed:
    raise SystemExit(
        "FAIL: QC structural assertions failed:\n  - "
        + "\n  - ".join(failed)
    )

match = re.search(
    r'cat > "fragment_pairs\.awk" <<\'AWK\'\n(.*?)\nAWK',
    bam_to_fragments,
    re.DOTALL,
)
if not match:
    raise SystemExit("FAIL: could not extract fragment_pairs.awk")
(Path.cwd() / ".fragment_pairs.test.awk").write_text(
    match.group(1).replace("\\$", "$") + "\n"
)
filtered_flag_match = re.search(
    r'cat > "filtered_flags\.awk" <<\'AWK\'\n(.*?)\nAWK',
    filtered_bam_qc,
    re.DOTALL,
)
if not filtered_flag_match:
    raise SystemExit("FAIL: could not extract filtered_flags.awk")
(Path.cwd() / ".filtered_flags.test.awk").write_text(
    filtered_flag_match.group(1).replace("\\$", "$") + "\n"
)
aggregators = re.findall(
    r"python <<'PY'\n(.*?)\nPY",
    multiqc,
    re.DOTALL,
)
aggregator = next(
    (source for source in aggregators if "duplicate peak-QC sample_id" in source),
    None,
)
if aggregator is None:
    raise SystemExit("FAIL: could not extract MultiQC aggregation program")
(Path.cwd() / ".multiqc_aggregate.test.py").write_text(aggregator + "\n")
library_formatter = re.search(
    r"process LIBRARY_QC_CUSTOM \{.*?"
    r"python - .*? <<'PY'\n(.*?)\nPY",
    multiqc,
    re.DOTALL,
)
if not library_formatter:
    raise SystemExit("FAIL: could not extract library-QC formatter")
(Path.cwd() / ".library_qc_formatter.test.py").write_text(
    library_formatter.group(1).replace("\\\\", "\\") + "\n"
)
motif_formatter = re.search(
    r"process MOTIF_QC_CUSTOM \{.*?"
    r"python - .*? <<'PY'\n(.*?)\nPY",
    multiqc,
    re.DOTALL,
)
if not motif_formatter:
    raise SystemExit("FAIL: could not extract motif-QC formatter")
(Path.cwd() / ".motif_qc_formatter.test.py").write_text(
    motif_formatter.group(1).replace("\\\\", "\\") + "\n"
)
gtf_prepare = re.search(
    r"def prepareTss = annotation_mode == 'gtf' \? '''\n(.*?)\n    ''' : '''",
    tss,
    re.DOTALL,
)
if not gtf_prepare:
    raise SystemExit("FAIL: could not extract strand-aware GTF conversion")
(Path.cwd() / ".tss_gtf.test.sh").write_text(
    "set -euo pipefail\n"
    + gtf_prepare.group(1).replace("__OUTPUT_STEM__", "strand_test")
    + "\n"
)
PY

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/nanocut-qc.XXXXXX")
trap 'rm -rf -- "${tmp_dir:?}" .fragment_pairs.test.awk .filtered_flags.test.awk .multiqc_aggregate.test.py .library_qc_formatter.test.py .motif_qc_formatter.test.py .tss_gtf.test.sh' EXIT
mkdir -p "$tmp_dir/direct"

cat > "$tmp_dir/direct/name_sorted.sam" <<'EOF'
pairA	99	chrMini	11	60	10M	=	41	40	AAAAAAAAAA	IIIIIIIIII
pairA	147	chrMini	41	60	10M	=	11	-40	TTTTTTTTTT	IIIIIIIIII
pairB	147	chrMini	81	60	5M2D5M	=	61	-30	CCCCCCCCCC	IIIIIIIIII
pairB	99	chrMini	61	60	10M	=	81	30	GGGGGGGGGG	IIIIIIIIII
EOF

LC_ALL=C awk -f .fragment_pairs.test.awk \
  "$tmp_dir/direct/name_sorted.sam" \
  | LC_ALL=C sort -t $'\t' \
      -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n -k7,7 \
  > "$tmp_dir/direct/fragments.bedpe"

[[ $(wc -l < "$tmp_dir/direct/fragments.bedpe") -eq 2 ]]
awk -F $'\t' 'NF != 7 || $7 == "" || $7 == "." { exit 1 }' \
  "$tmp_dir/direct/fragments.bedpe"
grep -F $'chrMini\t10\t20\tchrMini\t40\t50\tpairA' \
  "$tmp_dir/direct/fragments.bedpe" >/dev/null
grep -F $'chrMini\t60\t70\tchrMini\t80\t92\tpairB' \
  "$tmp_dir/direct/fragments.bedpe" >/dev/null

cat > "$tmp_dir/direct/two_first_mates.sam" <<'EOF'
badPair	99	chrMini	11	60	10M	=	41	40	AAAAAAAAAA	IIIIIIIIII
badPair	83	chrMini	41	60	10M	=	11	-40	TTTTTTTTTT	IIIIIIIIII
EOF
if LC_ALL=C awk -f .fragment_pairs.test.awk \
  "$tmp_dir/direct/two_first_mates.sam" \
  > "$tmp_dir/direct/bad.fragments.bedpe" 2> "$tmp_dir/bad-mates.log"; then
  printf 'FAIL: two 0x40 first mates unexpectedly formed a fragment\n' >&2
  exit 1
fi
grep -F 'expected one 0x40 and one 0x80 mate for badPair' \
  "$tmp_dir/bad-mates.log" >/dev/null

cat > "$tmp_dir/direct/both_and_neither.sam" <<'EOF'
xorPair	195	chrMini	11	60	10M	=	41	40	AAAAAAAAAA	IIIIIIIIII
xorPair	3	chrMini	41	60	10M	=	11	-40	TTTTTTTTTT	IIIIIIIIII
EOF
if LC_ALL=C awk -f .fragment_pairs.test.awk \
  "$tmp_dir/direct/both_and_neither.sam" \
  > "$tmp_dir/direct/xor.fragments.bedpe" 2> "$tmp_dir/bam-xor.log"; then
  printf 'FAIL: both-bit/neither-bit records unexpectedly formed a fragment\n' >&2
  exit 1
fi
grep -F 'each alignment must set exactly one of 0x40 and 0x80 for xorPair' \
  "$tmp_dir/bam-xor.log" >/dev/null

if LC_ALL=C awk -f .filtered_flags.test.awk \
  "$tmp_dir/direct/both_and_neither.sam" \
  > /dev/null 2> "$tmp_dir/filtered-xor.log"; then
  printf 'FAIL: filtered BAM validator accepted both-bit/neither-bit records\n' >&2
  exit 1
fi
grep -F 'each filtered alignment must set exactly one of 0x40 and 0x80' \
  "$tmp_dir/filtered-xor.log" >/dev/null

cat > "$tmp_dir/direct/final.broadPeak" <<'EOF'
chrMini	15	45	peak_1	50	.	4	2	1
chrMini	42	65	peak_2	60	.	5	3	2
EOF

python3 bin/peak_qc.py \
  --fragments "$tmp_dir/direct/fragments.bedpe" \
  --peaks "$tmp_dir/direct/final.broadPeak" \
  --sample-id TARGET \
  --output-prefix "$tmp_dir/direct/TARGET.peak_qc"

python3 - "$tmp_dir/direct/TARGET.peak_qc.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    metrics = json.load(handle)
assert metrics["sample_id"] == "TARGET"
assert metrics["total_fragments"] == 2
assert metrics["fragments_in_peaks"] == 2
assert metrics["frip"] == 1.0
assert metrics["peak_count"] == 2
PY

: > "$tmp_dir/direct/empty.broadPeak"
python3 bin/peak_qc.py \
  --fragments "$tmp_dir/direct/fragments.bedpe" \
  --peaks "$tmp_dir/direct/empty.broadPeak" \
  --sample-id EMPTY \
  --output-prefix "$tmp_dir/direct/EMPTY.peak_qc"

python3 - "$tmp_dir/direct/EMPTY.peak_qc.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    metrics = json.load(handle)
assert metrics["peak_count"] == 0
assert metrics["fragments_in_peaks"] == 0
assert metrics["frip"] == 0.0
PY

mkdir -p "$tmp_dir/library_formatter"
cat > "$tmp_dir/library_formatter/flagstat.txt" <<'EOF'
100 + 0 in total (QC-passed reads + QC-failed reads)
80 + 0 mapped (80.00% : N/A)
70 + 0 properly paired (70.00% : N/A)
99 + 0 duplicates
EOF
cat > "$tmp_dir/library_formatter/alignment.stats.txt" <<'EOF'
SN	raw total sequences:	100
SN	reads mapped:	80
SN	reads properly paired:	70
IS	100	2	2	0	0
IS	200	1	1	0	0
EOF
cat > "$tmp_dir/library_formatter/idxstats.tsv" <<'EOF'
chrMini	1000	75	0
chrM	100	5	0
*	0	0	20
EOF
cat > "$tmp_dir/library_formatter/insert_size.tsv" <<'EOF'
record_type	insert_size	total_pairs	inward_pairs	outward_pairs	other_pairs
IS	100	2	2	0	0
IS	200	1	1	0	0
EOF
cat > "$tmp_dir/library_formatter/duplicate_metrics.json" <<'EOF'
{"READ": 100.0, "DUPLICATE TOTAL": 20, "ESTIMATED LIBRARY SIZE": 40.0}
EOF
cat > "$tmp_dir/library_formatter/filtered_bam_qc.tsv" <<'EOF'
metric	value
filtered_reads	60
filtered_fragments	30
EOF
(
  cd "$tmp_dir/library_formatter"
  python3 "$repo_root/.library_qc_formatter.test.py" \
    TARGET TARGET.library_qc.tsv TARGET.insert_size_distribution.tsv
)
python3 - "$tmp_dir/library_formatter/TARGET.library_qc.tsv" <<'PY'
import csv
import math
import sys

with open(sys.argv[1], encoding="utf-8", newline="") as handle:
    row = next(csv.DictReader(handle, delimiter="\t"))
assert row["sample_id"] == "TARGET"
assert float(row["duplicate_percent"]) == 20.0
assert float(row["duplicate_total"]) == 20.0
assert float(row["mapq_filtered_fraction"]) == 0.6
assert int(float(row["mapq_filtered_reads"])) == 60
assert int(float(row["mapq_filtered_fragments"])) == 30
assert int(float(row["insert_size_total_pairs"])) == 3
assert math.isclose(float(row["insert_size_mean"]), 400 / 3)
assert float(row["insert_size_min"]) == 100
assert float(row["insert_size_max"]) == 200
PY
grep -F $'TARGET\t100\t2' \
  "$tmp_dir/library_formatter/TARGET.insert_size_distribution.tsv" >/dev/null
grep -F $'TARGET\t200\t1' \
  "$tmp_dir/library_formatter/TARGET.insert_size_distribution.tsv" >/dev/null

cat > "$tmp_dir/library_formatter/duplicate_metrics.json" <<'EOF'
{"READ": 100, "DUPLICATE TOTAL": 20}
EOF
(
  cd "$tmp_dir/library_formatter"
  python3 "$repo_root/.library_qc_formatter.test.py" \
    TARGET no-estimate.library_qc.tsv no-estimate.insert.tsv
)
python3 - "$tmp_dir/library_formatter/no-estimate.library_qc.tsv" <<'PY'
import csv
import sys

with open(sys.argv[1], encoding="utf-8", newline="") as handle:
    row = next(csv.DictReader(handle, delimiter="\t"))
assert row["estimated_library_size"] == ""
assert float(row["duplicate_total"]) == 20.0
assert float(row["duplicate_percent"]) == 20.0
PY

cat > "$tmp_dir/library_formatter/duplicate_metrics.json" <<'EOF'
{"READ": 100, "ESTIMATED LIBRARY SIZE": 40}
EOF
if (
  cd "$tmp_dir/library_formatter"
  python3 "$repo_root/.library_qc_formatter.test.py" \
    TARGET missing.tsv missing.insert.tsv \
    > "$tmp_dir/missing-markdup.log" 2>&1
); then
  printf 'FAIL: missing DUPLICATE TOTAL unexpectedly accepted\n' >&2
  exit 1
fi
grep -F 'missing required markdup field: DUPLICATE TOTAL' \
  "$tmp_dir/missing-markdup.log" >/dev/null

cat > "$tmp_dir/library_formatter/duplicate_metrics.json" <<'EOF'
{"EXAMINED": true, "DUPLICATE TOTAL": 20, "ESTIMATED LIBRARY SIZE": 40}
EOF
if (
  cd "$tmp_dir/library_formatter"
  python3 "$repo_root/.library_qc_formatter.test.py" \
    TARGET invalid.tsv invalid.insert.tsv \
    > "$tmp_dir/invalid-markdup.log" 2>&1
); then
  printf 'FAIL: boolean markdup count unexpectedly accepted\n' >&2
  exit 1
fi
grep -F 'markdup field EXAMINED must be a JSON number' \
  "$tmp_dir/invalid-markdup.log" >/dev/null

mkdir -p "$tmp_dir/motif_formatter"
cat > "$tmp_dir/motif_formatter/expected_motif_qc.tsv" <<'EOF'
metric	value
status	pass
best_motif_id	MA0139.1
best_adjusted_p_value	0.001
EOF
(
  cd "$tmp_dir/motif_formatter"
  python3 "$repo_root/.motif_qc_formatter.test.py" \
    TARGET TARGET.motif_qc.tsv
)
grep -F $'TARGET\tpass\tMA0139.1\t0.001' \
  "$tmp_dir/motif_formatter/TARGET.motif_qc.tsv" >/dev/null

mkdir -p \
  "$tmp_dir/aggregate/peak_qc01" \
  "$tmp_dir/aggregate/motif01" \
  "$tmp_dir/aggregate/insert01" \
  "$tmp_dir/aggregate/multiqc_custom_content"
cp "$tmp_dir/direct/TARGET.peak_qc.tsv" \
  "$tmp_dir/aggregate/peak_qc01/TARGET.peak_qc.tsv"
cp "$tmp_dir/motif_formatter/TARGET.motif_qc.tsv" \
  "$tmp_dir/aggregate/motif01/TARGET.motif_qc.tsv"
cp "$tmp_dir/library_formatter/TARGET.insert_size_distribution.tsv" \
  "$tmp_dir/aggregate/insert01/TARGET.insert_size_distribution.tsv"
(
  cd "$tmp_dir/aggregate"
  annotation_status=skipped_no_annotation \
    python3 "$repo_root/.multiqc_aggregate.test.py"
)
grep -F $'TARGET\t2\t' \
  "$tmp_dir/aggregate/combined_target_qc.tsv" >/dev/null
grep -F $'TARGET\tpass' \
  "$tmp_dir/aggregate/multiqc_custom_content/nanocut_motif_qc_mqc.tsv" \
  >/dev/null
grep -F $'insert_size\tTARGET' \
  "$tmp_dir/aggregate/multiqc_custom_content/nanocut_insert_size_mqc.tsv" \
  >/dev/null
grep -F $'100\t2' \
  "$tmp_dir/aggregate/multiqc_custom_content/nanocut_insert_size_mqc.tsv" \
  >/dev/null
grep -F $'<run>\tnone\tskipped_no_annotation' \
  "$tmp_dir/aggregate/multiqc_custom_content/nanocut_tss_qc_mqc.tsv" \
  >/dev/null

mkdir -p "$tmp_dir/aggregate/motif02"
cp "$tmp_dir/motif_formatter/TARGET.motif_qc.tsv" \
  "$tmp_dir/aggregate/motif02/TARGET.motif_qc.tsv"
if (
  cd "$tmp_dir/aggregate"
  annotation_status=skipped_no_annotation \
    python3 "$repo_root/.multiqc_aggregate.test.py" \
      > "$tmp_dir/duplicate-motif.log" 2>&1
); then
  printf 'FAIL: duplicate motif sample_id unexpectedly aggregated\n' >&2
  exit 1
fi
grep -F 'duplicate motif sample_id: TARGET' \
  "$tmp_dir/duplicate-motif.log" >/dev/null
rm -rf "$tmp_dir/aggregate/motif02"

mkdir -p "$tmp_dir/aggregate/peak_qc02"
cp "$tmp_dir/direct/TARGET.peak_qc.tsv" \
  "$tmp_dir/aggregate/peak_qc02/TARGET.peak_qc.tsv"
if (
  cd "$tmp_dir/aggregate"
  annotation_status=skipped_no_annotation \
    python3 "$repo_root/.multiqc_aggregate.test.py" \
      > "$tmp_dir/duplicate-summary.log" 2>&1
); then
  printf 'FAIL: duplicate peak-QC sample_id unexpectedly aggregated\n' >&2
  exit 1
fi
grep -F 'duplicate peak-QC sample_id: TARGET' \
  "$tmp_dir/duplicate-summary.log" >/dev/null

mkdir -p "$tmp_dir/tss/annotation"
cat > "$tmp_dir/tss/annotation/source.annotation" <<'EOF'
chrMini	source	transcript	11	20	.	+	.	gene_id "plus";
chrMini	source	transcript	31	50	.	-	.	gene_id "minus";
EOF
(
  cd "$tmp_dir/tss"
  bash "$repo_root/.tss_gtf.test.sh"
)
grep -F $'chrMini\t10\t11\ttss_1\t0\t+' \
  "$tmp_dir/tss/strand_test.tss.bed" >/dev/null
grep -F $'chrMini\t49\t50\ttss_2\t0\t-' \
  "$tmp_dir/tss/strand_test.tss.bed" >/dev/null

if ! command -v nextflow >/dev/null 2>&1; then
  printf '%s\n' \
    'SKIP: Nextflow runtime unavailable; strong static checks and direct mate-flag/FRiP/empty-peak/library-metrics/motif/TSS/summary checks passed, but no DSL2 execution was run.'
  exit 0
fi

mkdir -p "$tmp_dir/fakebin" "$tmp_dir/runtime/input"

cat > "$tmp_dir/fakebin/samtools" <<'PY'
#!/usr/bin/env python3
import pathlib
import sys

args = sys.argv[1:]
if args == ["--version"]:
    print("samtools 1.20")
elif args and args[0] == "sort":
    path = pathlib.Path(args[-1])
    sys.stdout.write(path.read_text())
elif args and args[0] == "view":
    if "-c" in args:
        records = [
            line.split("\t")
            for line in pathlib.Path(args[-1]).read_text().splitlines()
            if line and not line.startswith("@")
        ]
        required_flag = None
        if "-f" in args:
            required_flag = int(args[args.index("-f") + 1])
        count = sum(
            required_flag is None
            or (int(fields[1]) & required_flag) == required_flag
            for fields in records
        )
        print(count)
    else:
        if args[-1] == "-":
            sys.stdout.write(sys.stdin.read())
        else:
            sys.stdout.write(pathlib.Path(args[-1]).read_text())
else:
    raise SystemExit(f"unsupported fake samtools arguments: {args}")
PY
chmod +x "$tmp_dir/fakebin/samtools"

cat > "$tmp_dir/fakebin/multiqc" <<'PY'
#!/usr/bin/env python3
import pathlib
import sys

if sys.argv[1:] == ["--version"]:
    print("multiqc, version 1.25.2")
    raise SystemExit(0)
pathlib.Path("multiqc_report.html").write_text("<html>fake MultiQC</html>\n")
pathlib.Path("multiqc_data").mkdir(exist_ok=True)
pathlib.Path("multiqc_data/multiqc_sources.txt").write_text("fake\n")
PY
chmod +x "$tmp_dir/fakebin/multiqc"

cp "$tmp_dir/direct/name_sorted.sam" "$tmp_dir/runtime/input/IgG.bam"
cp "$tmp_dir/direct/name_sorted.sam" "$tmp_dir/runtime/input/TARGET.bam"
: > "$tmp_dir/runtime/input/IgG.bam.bai"
: > "$tmp_dir/runtime/input/TARGET.bam.bai"
cp "$tmp_dir/direct/final.broadPeak" "$tmp_dir/runtime/input/TARGET.broadPeak"
: > "$tmp_dir/runtime/input/IgG.bw"
: > "$tmp_dir/runtime/input/TARGET.bw"

for sample in IgG TARGET; do
  {
    printf '4 + 0 in total (QC-passed reads + QC-failed reads)\n'
    printf '4 + 0 mapped (100.00%% : N/A)\n'
    printf '4 + 0 properly paired (100.00%% : N/A)\n'
    printf '0 + 0 duplicates\n'
  } > "$tmp_dir/runtime/input/${sample}.flagstat.txt"
  {
    printf 'SN\traw total sequences:\t4\n'
    printf 'SN\treads mapped:\t4\n'
    printf 'SN\treads properly paired:\t4\n'
    printf 'IS\t40\t2\t2\t0\t0\n'
  } > "$tmp_dir/runtime/input/${sample}.stats.txt"
  printf 'chrMini\t100\t4\t0\nchrM\t20\t0\t0\n*\t0\t0\t0\n' \
    > "$tmp_dir/runtime/input/${sample}.idxstats.tsv"
  printf 'record_type\tinsert_size\ttotal_pairs\tinward_pairs\toutward_pairs\tother_pairs\n' \
    > "$tmp_dir/runtime/input/${sample}.insert_size.tsv"
  printf 'IS\t40\t2\t2\t0\t0\n' \
    >> "$tmp_dir/runtime/input/${sample}.insert_size.tsv"
  printf '{"READ": 4, "DUPLICATE TOTAL": 0, "ESTIMATED LIBRARY SIZE": 2}\n' \
    > "$tmp_dir/runtime/input/${sample}.duplicate_metrics.json"
done

cat > "$tmp_dir/runtime/input/demux.json" <<'EOF'
{"total_reads": 4, "assigned_reads": 4, "ambiguous_reads": 0, "unassigned_reads": 0, "assigned_fraction": 1.0, "ambiguous_fraction": 0.0, "unassigned_fraction": 0.0}
EOF
printf 'metric\tvalue\ntotal_reads\t4\nassigned_fraction\t1.0\n' \
  > "$tmp_dir/runtime/input/demux.tsv"
cat > "$tmp_dir/runtime/input/TARGET.motif_qc.tsv" <<'EOF'
metric	value
status	pass
best_motif_id	MA0139.1
best_adjusted_p_value	0.001
EOF

cat > "$tmp_dir/runtime/main.nf" <<EOF
nextflow.enable.dsl = 2

include { QC } from '${repo_root}/subworkflows/local/qc'

workflow {
    filtered_bams = Channel.of(
        tuple(
            [sample_id: 'IgG', is_control: true],
            file('${tmp_dir}/runtime/input/IgG.bam'),
            file('${tmp_dir}/runtime/input/IgG.bam.bai')
        ),
        tuple(
            [sample_id: 'TARGET', is_control: false],
            file('${tmp_dir}/runtime/input/TARGET.bam'),
            file('${tmp_dir}/runtime/input/TARGET.bam.bai')
        )
    )
    final_peaks = Channel.of(
        tuple(
            [sample_id: 'TARGET', is_control: false],
            file('${tmp_dir}/runtime/input/TARGET.broadPeak')
        )
    )
    coverage = Channel.of(
        tuple([sample_id: 'IgG', is_control: true], file('${tmp_dir}/runtime/input/IgG.bw')),
        tuple([sample_id: 'TARGET', is_control: false], file('${tmp_dir}/runtime/input/TARGET.bw'))
    )
    library_metrics = Channel.of(
        tuple(
            [sample_id: 'IgG', is_control: true],
            file('${tmp_dir}/runtime/input/IgG.flagstat.txt'),
            file('${tmp_dir}/runtime/input/IgG.stats.txt'),
            file('${tmp_dir}/runtime/input/IgG.idxstats.tsv'),
            file('${tmp_dir}/runtime/input/IgG.insert_size.tsv'),
            file('${tmp_dir}/runtime/input/IgG.duplicate_metrics.json')
        ),
        tuple(
            [sample_id: 'TARGET', is_control: false],
            file('${tmp_dir}/runtime/input/TARGET.flagstat.txt'),
            file('${tmp_dir}/runtime/input/TARGET.stats.txt'),
            file('${tmp_dir}/runtime/input/TARGET.idxstats.tsv'),
            file('${tmp_dir}/runtime/input/TARGET.insert_size.tsv'),
            file('${tmp_dir}/runtime/input/TARGET.duplicate_metrics.json')
        )
    )
    demux_metrics = Channel.of(
        tuple(
            [library_id: 'LIB'],
            file('${tmp_dir}/runtime/input/demux.json'),
            file('${tmp_dir}/runtime/input/demux.tsv')
        )
    )
    QC(
        filtered_bams,
        final_peaks,
        coverage,
        library_metrics,
        demux_metrics,
        Channel.empty(),
        Channel.of(
            tuple(
                [sample_id: 'TARGET', is_control: false],
                file('${tmp_dir}/runtime/input/TARGET.motif_qc.tsv')
            )
        ),
        Channel.empty(),
        Channel.empty()
    )
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
  PATH="${repo_root}/bin:$tmp_dir/fakebin:$PATH" \
    nextflow run main.nf -c nextflow.config -ansi-log false \
      -with-trace "$tmp_dir/runtime/trace.txt"
)

[[ -f "$tmp_dir/runtime/results/qc/peaks/TARGET/TARGET.peak_qc.json" ]]
[[ ! -e "$tmp_dir/runtime/results/qc/peaks/IgG/IgG.peak_qc.json" ]]
[[ -f "$tmp_dir/runtime/results/reports/summary/combined_target_qc.tsv" ]]
[[ -f "$tmp_dir/runtime/results/reports/multiqc/multiqc_report.html" ]]
[[ $(grep -c 'BAM_TO_FRAGMENTS' "$tmp_dir/runtime/trace.txt") -eq 1 ]]
[[ $(grep -c 'PEAK_QC' "$tmp_dir/runtime/trace.txt") -eq 1 ]]
[[ $(grep -c 'FILTERED_BAM_QC' "$tmp_dir/runtime/trace.txt") -eq 2 ]]
[[ $(grep -c 'LIBRARY_QC_CUSTOM' "$tmp_dir/runtime/trace.txt") -eq 2 ]]
[[ $(grep -c 'MOTIF_QC_CUSTOM' "$tmp_dir/runtime/trace.txt") -eq 1 ]]
[[ $(grep -c 'TSS_ENRICHMENT' "$tmp_dir/runtime/trace.txt" || true) -eq 0 ]]
grep -F $'TARGET\t2\t' \
  "$tmp_dir/runtime/results/reports/summary/combined_target_qc.tsv" >/dev/null
grep -F $'TARGET\tpass\tMA0139.1\t0.001' \
  "$tmp_dir/runtime/results/reports/multiqc/multiqc_custom_content/nanocut_motif_qc_mqc.tsv" \
  >/dev/null

printf '%s\n' \
  'PASS: filtered fragment QC, target-only FRiP, filtered/library metrics, metadata-keyed motif QC, empty optional annotation, and MultiQC aggregation executed successfully.'
