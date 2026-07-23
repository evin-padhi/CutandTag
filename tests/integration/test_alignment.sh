#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"

required_files=(
  modules/local/bowtie2_build.nf
  modules/local/bowtie2_align.nf
  modules/local/samtools_sort_index.nf
  modules/local/samtools_filter.nf
  modules/local/samtools_metrics.nf
  modules/local/bamcoverage.nf
  subworkflows/local/align_qc.nf
)

for required_file in "${required_files[@]}"; do
  if [[ ! -f "$required_file" ]]; then
    printf 'FAIL: missing %s\n' "$required_file" >&2
    exit 1
  fi
done

python3 - <<'PY'
from pathlib import Path

root = Path.cwd()
sources = {
    path.name: path.read_text()
    for path in [
        root / "modules/local/bowtie2_build.nf",
        root / "modules/local/bowtie2_align.nf",
        root / "modules/local/samtools_sort_index.nf",
        root / "modules/local/samtools_filter.nf",
        root / "modules/local/samtools_metrics.nf",
        root / "modules/local/bamcoverage.nf",
        root / "subworkflows/local/align_qc.nf",
    ]
}

build = sources["bowtie2_build.nf"]
align = sources["bowtie2_align.nf"]
sort_index = sources["samtools_sort_index.nf"]
filter_bam = sources["samtools_filter.nf"]
metrics = sources["samtools_metrics.nf"]
coverage = sources["bamcoverage.nf"]
subworkflow = sources["align_qc.nf"]

bowtie_container = (
    "quay.io/biocontainers/bowtie2:2.5.4--he20e202_2"
)
samtools_container = (
    "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
)
deeptools_container = (
    "quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0"
)

checks = {
    "BOWTIE2_BUILD process is declared":
        "process BOWTIE2_BUILD" in build,
    "reference process uses the project-qualified Bowtie2 environment":
        'conda "${projectDir}/envs/bowtie2.yml"' in build,
    "reference process uses the pinned valid Bowtie2 image":
        f"container '{bowtie_container}'" in build,
    "reference process accepts FASTA or supplied-index mode":
        "val(reference_mode)" in build
        and "'fasta'" in build
        and "'index'" in build,
    "reference inputs are staged below a fixed directory":
        "stageAs: 'reference_input/*'" in build,
    "FASTA mode invokes bowtie2-build against a quoted staged path":
        "bowtie2-build" in build and '"${reference_file}"' in build,
    "supplied small and large indexes are canonicalized safely":
        all(token in build for token in (
            "*.1.bt2)", "*.rev.2.bt2)", "*.1.bt2l)", "*.rev.2.bt2l)",
            "for index_suffix in 1 2 3 4 rev.1 rev.2",
            'bt2_index/reference.${index_suffix}.${index_extension}',
        )),
    "reference process emits canonical index files and a keyed version file":
        "tuple val(reference_meta)" in build
        and "bt2_index/reference*.bt2*" in build
        and "versions" in build,
    "BOWTIE2_ALIGN process is declared":
        "process BOWTIE2_ALIGN" in align,
    "alignment uses the project-qualified Bowtie2 environment and image":
        'conda "${projectDir}/envs/bowtie2.yml"' in align
        and f"container '{bowtie_container}'" in align,
    "alignment stages reads and indexes with collision-proof aliases":
        all(token in align for token in (
            "path(r1, stageAs: 'reads_R1.fastq.gz')",
            "path(r2, stageAs: 'reads_R2.fastq.gz')",
            "path(index_files, stageAs: 'bt2_index/*')",
        )),
    "alignment uses only fixed quoted shell paths":
        '-x "bt2_index/reference"' in align
        and '-1 "reads_R1.fastq.gz"' in align
        and '-2 "reads_R2.fastq.gz"' in align
        and '-S "alignment.sam"' in align,
    "alignment captures the Bowtie2 summary separately":
        '2> "bowtie2.summary.txt"' in align,
    "alignment emits metadata-keyed SAM, summary, and version outputs":
        align.count("tuple val(meta)") >= 3
        and "bowtie2.summary.txt" in align
        and "versions" in align,
    "sample metadata is never interpolated into module shell scripts":
        all(
            "${meta." not in source.split("script:", 1)[1]
            for source in (align, sort_index, filter_bam, metrics, coverage)
        ),
    "SAMTOOLS_SORT_INDEX process is declared":
        "process SAMTOOLS_SORT_INDEX" in sort_index,
    "sort/index uses the project-qualified SAMtools environment and image":
        'conda "${projectDir}/envs/samtools.yml"' in sort_index
        and f"container '{samtools_container}'" in sort_index,
    "sort/index preserves every alignment in the primary analysis BAM":
        'samtools sort' in sort_index
        and "samtools view" not in sort_index
        and "analysis.bam" in sort_index,
    "sort/index emits metadata-keyed analysis BAM and BAI":
        "tuple val(meta)" in sort_index
        and 'path("analysis.bam")'
        in sort_index
        and 'path("analysis.bam.bai")' in sort_index,
    "SAMTOOLS_FILTER process is declared":
        "process SAMTOOLS_FILTER" in filter_bam,
    "filtering uses the project-qualified SAMtools environment and image":
        'conda "${projectDir}/envs/samtools.yml"' in filter_bam
        and f"container '{samtools_container}'" in filter_bam,
    "filtering validates MAPQ as a non-negative integer before interpolation":
        "Integer.parseInt" in filter_bam
        and filter_bam.index("Integer.parseInt") < filter_bam.index("samtools view"),
    "filtering requires proper pairs and excludes secondary/supplementary reads":
        "-f 2" in filter_bam and "-F 2304" in filter_bam,
    "filtering keeps duplicates for explicit downstream duplicate handling":
        "-F 3328" not in filter_bam and "-F 1280" not in filter_bam,
    "filtering indexes and emits the metadata-keyed filtered BAM":
        "samtools index" in filter_bam
        and "filtered.bam.bai" in filter_bam
        and "tuple val(meta)" in filter_bam,
    "SAMTOOLS_METRICS process is declared":
        "process SAMTOOLS_METRICS" in metrics,
    "metrics uses the project-qualified SAMtools environment and image":
        'conda "${projectDir}/envs/samtools.yml"' in metrics
        and f"container '{samtools_container}'" in metrics,
    "metrics emits flagstat, stats, idxstats, and insert-size rows":
        all(token in metrics for token in (
            "samtools flagstat",
            "samtools stats",
            "samtools idxstats",
            "insert_size.tsv",
            "$1 == \"IS\"",
        )),
    "metrics computes SAMtools-native duplicate and library-complexity stats":
        all(token in metrics for token in (
            "samtools collate",
            "samtools fixmate -m",
            "samtools markdup",
            "--json",
            "duplicate_metrics.json",
        )),
    "metrics emits metadata-keyed metric and version outputs":
        metrics.count("tuple val(meta)") >= 2
        and "versions" in metrics,
    "BAMCOVERAGE process is declared":
        "process BAMCOVERAGE" in coverage,
    "bamCoverage uses the project-qualified deepTools environment and image":
        'conda "${projectDir}/envs/deeptools.yml"' in coverage
        and f"container '{deeptools_container}'" in coverage,
    "bamCoverage uses the exact NanoScope coverage settings":
        all(token in coverage for token in (
            "--minMappingQuality",
            "--binSize 50",
            "--centerReads",
            "--smoothLength 250",
            "--normalizeUsing RPKM",
            "--ignoreDuplicates",
            "--extendReads",
        )),
    "bamCoverage receives the validated MAPQ and fixed staged BAM path":
        '"${mapq}"' in coverage
        and '--bam "coverage.filtered.bam"' in coverage,
    "bamCoverage emits metadata-keyed bigWig and versions":
        "tuple val(meta)" in coverage
        and "coverage.RPKM.bw" in coverage
        and "versions" in coverage,
    "alignment subworkflow includes every Task 7 process":
        all(name in subworkflow for name in (
            "BOWTIE2_BUILD",
            "BOWTIE2_ALIGN",
            "SAMTOOLS_SORT_INDEX",
            "SAMTOOLS_FILTER",
            "SAMTOOLS_METRICS",
            "BAMCOVERAGE",
        )),
    "subworkflow validates sample metadata as a safe filename contract":
        "SAFE_ID" in subworkflow
        and "sample_id" in subworkflow
        and "matches(SAFE_ID)" in subworkflow,
    "subworkflow resolves supplied small or large index files":
        ".bt2{,l}" in subworkflow
        and "bowtie2_index" in subworkflow,
    "subworkflow falls back to FASTA index construction":
        "fasta" in subworkflow and "reference_mode" in subworkflow,
    "subworkflow fans the canonical index out to every sample":
        ".combine(BOWTIE2_BUILD.out.index)" in subworkflow,
    "subworkflow preserves analysis BAM for metrics and downstream peak calls":
        "SAMTOOLS_METRICS(SAMTOOLS_SORT_INDEX.out.bam)" in subworkflow,
    "subworkflow uses only filtered BAM for coverage":
        "BAMCOVERAGE(SAMTOOLS_FILTER.out.bam" in subworkflow,
    "subworkflow exposes metadata-keyed BAMs, metrics, coverage, summaries, and versions":
        all(token in subworkflow for token in (
            "analysis_bam =",
            "filtered_bam =",
            "metrics =",
            "coverage =",
            "alignment_summary =",
            "versions =",
        )),
}

failed = [description for description, passed in checks.items() if not passed]
if failed:
    raise SystemExit(
        "FAIL: alignment structural assertions failed:\n  - "
        + "\n  - ".join(failed)
    )
PY

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/nanocut-alignment.XXXXXX")
trap 'rm -rf -- "${tmp_dir:?}"' EXIT

mkdir -p "$tmp_dir/supplied"
for suffix in 1 2 3 4 rev.1 rev.2; do
  : > "$tmp_dir/supplied/mini.${suffix}.bt2"
done

cat > "$tmp_dir/reference.fa" <<'EOF'
>chrMini
TTTTACGTCAGTACGATCGATGCTAGCAACCGGTTAACCGGTTAGCTAGCATCGATCGTACGATGAAAA
EOF

cat > "$tmp_dir/reads_R1.fastq" <<'EOF'
@pair-1/1
ACGTCAGTACGATCGA
+
IIIIIIIIIIIIIIII
@pair-2/1
AACCGGTTAACCGGTT
+
IIIIIIIIIIIIIIII
EOF

cat > "$tmp_dir/reads_R2.fastq" <<'EOF'
@pair-1/2
GCTAGCATCGATCGTA
+
IIIIIIIIIIIIIIII
@pair-2/2
ATCGATCGTACGATGA
+
IIIIIIIIIIIIIIII
EOF

python3 - "$tmp_dir/supplied/mini" "$tmp_dir/reference.fa" "$tmp_dir/reads_R1.fastq" "$tmp_dir/reads_R2.fastq" <<'PY'
import sys
from pathlib import Path

prefix = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
r1_path = Path(sys.argv[3])
r2_path = Path(sys.argv[4])
suffixes = ("1", "2", "3", "4", "rev.1", "rev.2")
small = [Path(f"{prefix}.{suffix}.bt2") for suffix in suffixes]
large = [Path(f"{prefix}.{suffix}.bt2l") for suffix in suffixes]
assert all(path.is_file() for path in small)
assert not any(path.exists() for path in large)
assert {path.name for path in small} == {
    "mini.1.bt2",
    "mini.2.bt2",
    "mini.3.bt2",
    "mini.4.bt2",
    "mini.rev.1.bt2",
    "mini.rev.2.bt2",
}

reference = "".join(
    line.strip()
    for line in reference_path.read_text().splitlines()
    if not line.startswith(">")
)

def fastq_sequences(path):
    lines = path.read_text().splitlines()
    assert len(lines) % 4 == 0
    return [lines[index + 1] for index in range(0, len(lines), 4)]

complement = str.maketrans("ACGT", "TGCA")
for sequence in fastq_sequences(r1_path) + fastq_sequences(r2_path):
    reverse_complement = sequence.translate(complement)[::-1]
    assert sequence in reference or reverse_complement in reference, sequence
PY

if ! command -v nextflow >/dev/null 2>&1; then
  printf 'SKIP: Nextflow runtime unavailable; strong static module/subworkflow checks and direct supplied-index fixture checks passed, but no Nextflow execution was run.\n'
  exit 0
fi

fake_bin="$tmp_dir/bin"
mkdir -p "$fake_bin"

cat > "$fake_bin/bowtie2-build" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --version ]]; then
  printf 'bowtie2-build version 2.5.4\n'
  exit 0
fi
while [[ ${1:-} == --threads ]]; do
  shift 2
done
reference=${1:?missing reference}
prefix=${2:?missing index prefix}
[[ -s $reference ]]
for suffix in 1 2 3 4 rev.1 rev.2; do
  printf 'synthetic index\n' > "${prefix}.${suffix}.bt2"
done
EOF

cat > "$fake_bin/bowtie2" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --version ]]; then
  printf '/synthetic/bowtie2-align-s version 2.5.4\n'
  exit 0
fi
output=
while (( $# > 0 )); do
  case "$1" in
    --threads|-x|-1|-2)
      shift 2
      ;;
    -S)
      output=$2
      shift 2
      ;;
    *)
      printf 'unexpected fake bowtie2 argument: %s\n' "$1" >&2
      exit 1
      ;;
  esac
done
cat > "$output" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrMini	LN:76
pair-1	99	chrMini	5	42	16M	=	45	56	ACGTCAGTACGATCGA	IIIIIIIIIIIIIIII
pair-1	147	chrMini	45	42	16M	=	5	-56	GCTAGCATCGATCGTA	IIIIIIIIIIIIIIII
SAM
printf '1 reads; of these:\n  1 (100.00%%) aligned concordantly exactly 1 time\n100.00%% overall alignment rate\n' >&2
EOF

cat > "$fake_bin/samtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --version ]]; then
  printf 'samtools 1.20\n'
  exit 0
fi
command=${1:?missing samtools command}
shift
case "$command" in
  sort)
    output=
    input=
    while (( $# > 0 )); do
      case "$1" in
        -@)
          shift 2
          ;;
        -u)
          shift
          ;;
        -o)
          output=$2
          shift 2
          ;;
        -)
          input=-
          shift
          ;;
        *)
          input=$1
          shift
          ;;
      esac
    done
    if [[ $input == - ]]; then
      cat > "$output"
    else
      cp -- "$input" "$output"
    fi
    ;;
  quickcheck)
    exit 0
    ;;
  index)
    args=()
    while (( $# > 0 )); do
      if [[ $1 == -@ ]]; then
        shift 2
      else
        args+=("$1")
        shift
      fi
    done
    : > "${args[1]}"
    ;;
  view)
    output=
    input=
    while (( $# > 0 )); do
      case "$1" in
        -@|-f|-F|-q)
          shift 2
          ;;
        -b)
          shift
          ;;
        -o)
          output=$2
          shift 2
          ;;
        *)
          input=$1
          shift
          ;;
      esac
    done
    cp -- "$input" "$output"
    ;;
  flagstat)
    printf '2 + 0 in total (QC-passed reads + QC-failed reads)\n2 + 0 properly paired (100.00%% : N/A)\n'
    ;;
  stats)
    printf 'SN\traw total sequences:\t2\nIS\t56\t1\t1\t0\t0\n'
    ;;
  idxstats)
    printf 'chrMini\t76\t2\t0\n*\t0\t0\t0\n'
    ;;
  collate)
    input=${!#}
    cat "$input"
    ;;
  fixmate)
    cat
    ;;
  markdup)
    stats=
    while (( $# > 0 )); do
      case "$1" in
        -f)
          stats=$2
          shift 2
          ;;
        -@|-O)
          shift 2
          ;;
        --json|-s)
          shift
          ;;
        *)
          shift
          ;;
      esac
    done
    printf '{"READ":2,"DUPLICATE TOTAL":0,"ESTIMATED LIBRARY SIZE":1}\n' > "$stats"
    ;;
  *)
    printf 'unexpected fake samtools command: %s\n' "$command" >&2
    exit 1
    ;;
esac
EOF

cat > "$fake_bin/bamCoverage" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --version ]]; then
  printf 'bamCoverage 3.5.5\n'
  exit 0
fi
output=
while (( $# > 0 )); do
  case "$1" in
    --bam|--numberOfProcessors|--minMappingQuality|--binSize|--smoothLength|--normalizeUsing)
      shift 2
      ;;
    --outFileName)
      output=$2
      shift 2
      ;;
    --centerReads|--ignoreDuplicates|--extendReads)
      shift
      ;;
    *)
      printf 'unexpected fake bamCoverage argument: %s\n' "$1" >&2
      exit 1
      ;;
  esac
done
printf 'synthetic bigWig\n' > "$output"
EOF
chmod +x "$fake_bin/bowtie2-build" "$fake_bin/bowtie2" "$fake_bin/samtools" "$fake_bin/bamCoverage"

harness="$repo_root/.nextflow-task7-test-$$.nf"
runtime_config="$tmp_dir/nextflow.config"
trace_file="$tmp_dir/trace.tsv"
trap 'rm -f -- "${harness:?}"; rm -rf -- "${tmp_dir:?}"' EXIT

cat > "$harness" <<'EOF'
nextflow.enable.dsl = 2

include { ALIGN_QC } from './subworkflows/local/align_qc'

workflow {
    reads_ch = Channel.of(
        tuple(
            [sample_id: 'sampleA', library_id: 'libA'],
            file(params.r1, checkIfExists: true),
            file(params.r2, checkIfExists: true)
        )
    )
    ALIGN_QC(
        reads_ch,
        params.fasta,
        params.bowtie2_index,
        params.min_mapq
    )
}
EOF

cat > "$runtime_config" <<EOF
params.outdir = '$tmp_dir/results'
params.bowtie2_index = ''
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
  --fasta "$tmp_dir/reference.fa" \
  --r1 "$tmp_dir/reads_R1.fastq" \
  --r2 "$tmp_dir/reads_R2.fastq" \
  --min_mapq 5 \
  -with-trace "$trace_file" \
  -ansi-log false

for process_name in BOWTIE2_BUILD BOWTIE2_ALIGN SAMTOOLS_SORT_INDEX SAMTOOLS_FILTER SAMTOOLS_METRICS BAMCOVERAGE; do
  observed=$(awk -F '\t' -v name="$process_name" \
    'NR > 1 && $3 ~ name { count++ } END { print count + 0 }' \
    "$trace_file")
  [[ $observed -eq 1 ]] || {
    printf 'FAIL: expected one %s task, observed %s\n' "$process_name" "$observed" >&2
    exit 1
  }
done

python3 - "$tmp_dir/results" <<'PY'
import json
import sys
from pathlib import Path

results = Path(sys.argv[1])
expected = [
    results / "alignment/sampleA/bowtie2.summary.txt",
    results / "alignment/sampleA/analysis.bam",
    results / "alignment/sampleA/analysis.bam.bai",
    results / "alignment/sampleA/filtered.bam",
    results / "alignment/sampleA/filtered.bam.bai",
    results / "coverage/sampleA/coverage.RPKM.bw",
    results / "qc/library/sampleA/flagstat.txt",
    results / "qc/library/sampleA/alignment.stats.txt",
    results / "qc/library/sampleA/idxstats.tsv",
    results / "qc/library/sampleA/insert_size.tsv",
    results / "qc/library/sampleA/duplicate_metrics.json",
]
for path in expected:
    assert path.is_file(), path

summary = (results / "alignment/sampleA/bowtie2.summary.txt").read_text()
assert "100.00% overall alignment rate" in summary
duplicates = json.loads(
    (results / "qc/library/sampleA/duplicate_metrics.json").read_text()
)
assert duplicates["ESTIMATED LIBRARY SIZE"] == 1
assert (
    results / "pipeline_info/bowtie2_index/reference.rev.2.bt2"
).is_file()
PY

printf 'PASS: miniature reference alignment emitted indexed analysis/filtered BAMs, complete SAMtools metrics, and NanoScope coverage.\n'
