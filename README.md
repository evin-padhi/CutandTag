# Bulk nano-CUT&Tag Nextflow pipeline

This Nextflow DSL2 workflow demultiplexes paired-end bulk nano-CUT&Tag reads
with an 8-base I2 barcode, aligns each derived sample, calls matched-IgG
NanoScope-compatible broad peaks, calculates library/peak QC, and optionally
tests expected and de novo motifs. The repository includes the exact
six-library manifest mapping in `assets/samples.example.csv`.

## Prerequisites

- A POSIX shell and Java suitable for the installed Nextflow release.
- Nextflow 23.10.0 or later (enforced by `manifest.nextflowVersion`).
- One supported dependency backend:
  - Conda, Mamba, or Micromamba with `-profile conda`; or
  - Docker with `-profile docker`.
- A reference FASTA or complete Bowtie2 index.
- The input manifest and readable R1, R2, and I2 FASTQs.
- A MACS2 genome-size shortcut such as `hs`, or a positive effective genome
  size integer.

The `conda` profile is the most direct portable installation path. Nextflow
automatically creates each process environment from the repository-local,
exactly pinned YAML file when it is first needed; users should not create a
single combined environment by hand. The Docker profile uses the
process-specific pinned images declared in the modules.

Nextflow runtime execution was not verified in this workspace. Unit tests,
static assertions, and direct synthetic fixtures were verified, but Nextflow
was unavailable; therefore DSL2 execution, Conda solving, container pulls,
real bioinformatics-tool compatibility, and cache/resume behavior still need
verification on a supported runtime.

## Dependency environments and caches

The canonical environments are:

| File | Direct pinned tool | Used for |
|---|---|---|
| `envs/python.yml` | Python 3.12.3 | Manifest validation, demultiplexing, custom QC, pipeline metadata |
| `envs/fastqc.yml` | FastQC 0.12.1 | Read QC |
| `envs/bowtie2.yml` | Bowtie2 2.5.4 | Index construction and alignment |
| `envs/samtools.yml` | SAMtools 1.20 | BAM sorting/filtering/indexing, metrics, fragments |
| `envs/deeptools.yml` | deepTools 3.5.5 | RPKM bigWig and optional TSS enrichment |
| `envs/macs2.yml` | MACS2 2.2.9.1 | Broad primary and narrow motif-only peak calls |
| `envs/bedtools.yml` | BEDTools 2.31.1 | Blacklist filtering and motif sequence/background intervals |
| `envs/meme.yml` | MEME Suite 5.5.7 | AME, STREME, and FIMO |
| `envs/multiqc.yml` | MultiQC 1.25.2 | Consolidated report and summary tables |

Each environment uses `conda-forge` before `bioconda`. On a local workstation,
the automatic cache defaults to `.nextflow-conda` under the launch directory:

```bash
nextflow run main.nf \
  -profile conda \
  --input assets/samples.example.csv \
  --fasta /references/GRCh38.fa \
  --macs_genome_size hs \
  --outdir results
```

For a cluster, place both the Conda cache and work directory on storage visible
with the same path from every compute node:

```bash
export NXF_CONDA_CACHEDIR=/shared/nanocut/nextflow-conda
export NXF_HOME=/shared/nanocut/nextflow-home

nextflow run main.nf \
  -profile conda \
  -work-dir /shared/nanocut/work \
  --input /shared/project/samples.csv \
  --fasta /shared/references/GRCh38.fa \
  --macs_genome_size hs \
  --outdir /shared/project/results
```

`NXF_CONDA_CACHEDIR` overrides `conda.cacheDir`; without it the local automatic
cache is used. `NXF_HOME` is optional and controls Nextflow's own home/cache
location, not the pipeline result directory. The shared directories must be
writable during initial environment creation and readable by all jobs.
Mamba/Micromamba selection may be supplied in a user Nextflow configuration;
the repository YAML files remain the dependency source of truth.

The offline `test` profile modifies `PATH`, `LC_ALL`, `LANG`, and
`PYTHONDONTWRITEBYTECODE` only inside its processes so repository fake tools
are deterministic. Those are test-only environment variables and are not
production configuration knobs.

## Manifest

The CSV has one row per barcode-derived sample and this exact header:

```csv
sample_id,library_id,input_group,barcode,assay_target,is_control,control_id,expected_motif,r1,r2,i2
```

The columns mean:

- `sample_id`: unique derived-sample identifier. It must begin with a letter or
  digit and contain only letters, digits, `.`, `_`, or `-`.
- `library_id`: physical-library identifier, with the same character rule.
  Rows from one library repeat the same FASTQ paths intentionally; the
  physical library is streamed only once.
- `input_group`: biological input group used to require same-input matched
  controls, for example `25K`, `50K`, or `100K`.
- `barcode`: expected I2 sequence in FASTQ orientation. All manifest barcodes
  must be non-empty and equal length, and barcodes within a library must be
  unique.
- `assay_target`: one of `IgG`, `CTCF`, `GATA1`, or `RUNX1`.
- `is_control`: `true`, `1`, or `yes` for controls; `false`, `0`, or `no` for
  targets, case-insensitively.
- `control_id`: matched IgG `sample_id`. It is blank for controls and required
  for targets. The referenced row must be an IgG control with the same
  `input_group`; it may come from another physical library.
- `expected_motif`: regular expression matched against motif ID and alternate
  name. It is blank for controls and required for every target.
- `r1`, `r2`, `i2`: synchronized FASTQ paths. Plain or gzip-compressed content
  is accepted. Relative paths resolve from the manifest directory, and all
  three paths must be identical across rows sharing a `library_id`.

Validation stops before expensive work for missing/blank required fields,
unsupported targets, invalid identifiers/booleans, missing FASTQs, duplicate
sample IDs, duplicate within-library barcodes, unequal barcode lengths,
inconsistent library FASTQ paths, invalid controls, or cross-input control
links. R1/R2/I2 record counts and normalized read identifiers are also checked
during streaming demultiplexing. A unique closest barcode within
`--barcode_mismatches` is assigned; ties are ambiguous and out-of-threshold
reads are unassigned.

### Exact six-library mapping

| Physical library | Input | `TATAGCCT` (barcode A) | `ATAGAGGC` (barcode B) | Matched IgG for targets |
|---|---:|---|---|---|
| NX702 | 25K | NX702_IgG | NX702_CTCF | NX702_IgG |
| NX703 | 50K | NX703_IgG | NX703_CTCF | NX703_IgG |
| NX701 | 100K | NX701_IgG | NX701_CTCF | NX701_IgG |
| NX704 | 25K | NX704_GATA1 | NX704_RUNX1 | NX702_IgG |
| NX705 | 50K | NX705_GATA1 | NX705_RUNX1 | NX703_IgG |
| NX706 | 100K | NX706_GATA1 | NX706_RUNX1 | NX701_IgG |

Thus NX701–NX703 use barcode A for IgG and barcode B for CTCF, while
NX704–NX706 use barcode A for GATA1 and barcode B for RUNX1. Copy
`assets/samples.example.csv`, update its FASTQ paths, and retain the listed
control relationships. The expected motif expressions are `CTCF`, `GATA1`,
and `RUNX1` for their respective targets.

## Reference, index, and optional annotations

Use one of:

- `--fasta /path/genome.fa`: a readable FASTA whose first nonblank line is a
  `>` header. The pipeline builds and publishes a Bowtie2 index.
- `--bowtie2_index /path/index/prefix`: a prefix for exactly one complete
  six-file `.bt2` or `.bt2l` set:
  `1`, `2`, `3`, `4`, `rev.1`, and `rev.2`. Supply the prefix only, without a
  suffix. A supplied index takes precedence for alignment.

`--fasta` is still required when `--motif_db` is supplied because motif
windows must be extracted from the genome. Optional `--blacklist` is a BED
file removed from final broad peaks and motif inputs. Optional `--tss_bed`
must be strand-aware BED6 and takes precedence over `--gtf`; otherwise
transcript features in the GTF are converted to strand-aware single-base TSS
records.

## Run the pipeline

Conda:

```bash
nextflow run main.nf \
  -profile conda \
  --input samples.csv \
  --fasta GRCh38.fa \
  --outdir results \
  --blacklist hg38-blacklist.v2.bed \
  --gtf gencode.annotation.gtf \
  --motif_db JASPAR2026_CORE_vertebrates.meme \
  --macs_genome_size hs
```

Docker:

```bash
nextflow run main.nf \
  -profile docker \
  --input samples.csv \
  --bowtie2_index /references/grch38/bowtie2/grch38 \
  --outdir results \
  --macs_genome_size 2913022398
```

The second example omits motif analysis; add both `--fasta` and `--motif_db`
to enable it. A user config can override executor/resources without modifying
the pipeline:

```bash
nextflow run main.nf \
  -profile conda \
  -c cluster.config \
  --input samples.csv \
  --fasta GRCh38.fa \
  --macs_genome_size hs \
  --outdir results
```

The deterministic offline orchestration fixture is:

```bash
nextflow run main.nf -profile test
```

It uses tiny bundled data and fake command implementations; it is not a
biological validation and was not run here because Nextflow was unavailable.

## Parameters

Every pipeline CLI parameter is listed below. Nextflow runtime options such as
`-profile`, `-c`, `-work-dir`, and `-resume` use one leading dash; pipeline
parameters use two.

| Parameter | Required/default | Meaning |
|---|---|---|
| `--input` | Required | CSV manifest, one row per barcode-derived sample. |
| `--fasta` | Required unless `--bowtie2_index` is supplied | Reference FASTA; also required for motif analysis. |
| `--bowtie2_index` | Required unless `--fasta` is supplied | Complete Bowtie2 index prefix; takes precedence for alignment. |
| `--outdir` | `results` | Result directory. |
| `--blacklist` | Absent | Optional BED regions removed from final broad peaks and motif inputs. |
| `--gtf` | Absent | Optional GTF for transcript-derived, strand-aware TSS positions. |
| `--tss_bed` | Absent | Optional BED6 TSS file; takes precedence over `--gtf`. |
| `--motif_db` | Absent | MEME-format database beginning with `MEME version` and containing at least one `MOTIF` record; enables motif analysis. |
| `--barcode_mismatches` | `0` | Non-negative I2 Hamming-distance threshold, smaller than barcode length. |
| `--allow_empty` | `false` | Permit a derived sample with zero assigned read pairs for diagnostic runs. |
| `--min_mapq` | `5` | Integer 0–255 used for the filtered BAM, coverage, and fragment QC. |
| `--macs_genome_size` | Required | Positive effective genome size or path-safe MACS2 shortcut such as `hs`. |
| `--macs_llocal` | Fixed `100000` | NanoScope-compatible MACS2 local lambda window; other values are rejected. |
| `--macs_keep_dup` | Fixed `1` | NanoScope-compatible MACS2 duplicate setting; other values are rejected. |
| `--macs_broad_cutoff` | Fixed `0.1` | Primary broad-peak cutoff; other values are rejected. |
| `--macs_max_gap` | Fixed `1000` | Primary broad-peak maximum gap; other values are rejected. |
| `--motif_use_narrow_peaks` | `true` | Use matched-control narrow summits for motif windows. If false, use final broad-peak midpoints. |
| `--motif_window` | `200` | Positive total motif-window width in bp. Windows shift at reference edges; a shorter contig yields its full length. |

Boolean CLI values must reach Nextflow as booleans (`true` or `false`).
Optional features are skipped only when their path parameter is absent; a
supplied invalid path fails validation.

## Outputs

The principal layout under `--outdir` is:

```text
results/
  demultiplex/<library_id>/
  fastqc/<sample_id>/
  alignment/<sample_id>/
  coverage/<sample_id>/
  peaks/<sample_id>/broad/raw/
  peaks/<sample_id>/broad/final/
  peaks/<sample_id>/narrow_motif_qc/
  qc/library/<sample_id>/
  qc/fragments/<sample_id>/
  qc/peaks/<sample_id>/
  qc/tss/<sample_id>/
  motifs/<sample_id>/
  reports/multiqc/
  reports/summary/
  pipeline_info/
```

- `demultiplex/` contains per-sample paired FASTQs plus library JSON/TSV
  assignment metrics.
- `fastqc/` contains paired FastQC HTML/ZIP outputs.
- `alignment/` contains the primary analysis BAM/BAI, MAPQ-filtered BAM/BAI,
  Bowtie2 summary, and process version records.
- `coverage/` contains RPKM `coverage.RPKM.bw` files generated with MAPQ
  filtering, 50-bp bins, centered/extended reads, 250-bp smoothing, and
  duplicate ignoring.
- `peaks/` retains the matched-IgG raw broadPeak/gappedPeak/XLS/log, the
  blacklist-filtered `final.broadPeak`, and optional motif-only narrowPeak and
  summits.
- `qc/library/` contains SAMtools flagstat/stats/idxstats, insert sizes,
  markdup metrics, filtered-read/fragment counts, and custom summary inputs.
- `qc/fragments/` and `qc/peaks/` contain target-only BEDPE fragments,
  fragment-based FRiP, peak metrics, width histograms, and per-peak counts.
- `qc/tss/` contains optional TSS BED, matrix, table, profile, and status.
- `motifs/` contains foreground/background intervals and FASTAs, AME known
  enrichment, STREME de novo discovery, FIMO scans, status/log files, and
  `expected_motif_qc` JSON/TSV/position outputs.
- `reports/multiqc/multiqc_report.html` is the top-level report.
  `reports/summary/combined_target_qc.tsv` is the combined target broad-peak
  summary.
- `pipeline_info/` contains the normalized manifest, validated parameter JSON,
  software versions, completion summary, built index when applicable,
  execution report, timeline, trace, and DAG.

Empty target peak sets are valid outputs with zero/NA QC and recorded motif
skip status. IgG libraries receive read/alignment/library QC but are not
peak-called against themselves and do not receive target FRiP by default.

## QC interpretation

- **Demultiplexing:** `assigned_fraction`, `ambiguous_fraction`, and
  `unassigned_fraction` use all synchronized input read pairs as denominator.
  Ambiguous means two or more expected barcodes tied at the closest permitted
  Hamming distance. Unassigned means none passed the threshold.
- **Alignment and pairing:** `mapped_percent` and
  `properly_paired_percent` are the percentages reported by `samtools
  flagstat` on the primary analysis BAM.
- **MAPQ-filtered fraction:** the filtered BAM retains properly paired primary
  alignments (`samtools view -f 2 -F 2304`) at or above `--min_mapq`.
  `mapq_filtered_fraction` is filtered alignment records divided by SAMtools
  `raw total sequences`; two validated mate records constitute one filtered
  fragment.
- **Duplicate fraction:** `duplicate_percent` is SAMtools markdup
  `DUPLICATE TOTAL` divided by `EXAMINED` (or `READ` for compatible output),
  reported as a percentage. `estimated_library_size` is retained. Duplicates
  are measured, not removed from the primary analysis BAM. MACS2 intentionally
  uses `--keep-dup 1`; coverage independently uses `--ignoreDuplicates`.
- **Mitochondrial fraction:** mapped records on `chrM`, `MT`, or `M` divided
  by mapped records on all named reference sequences, reported as a
  percentage.
- **FRiP:** for each non-control sample, the numerator is the number of unique,
  properly paired fragments from the filtered BAM overlapping at least one
  final (blacklist-filtered when applicable) broad peak. The denominator is
  all unique properly paired fragments in that filtered BAM. A fragment is
  counted once even if both mates or multiple peaks overlap.
- **Peak QC:** includes count, union-covered bases, width min/mean/median/max
  and quartiles, width histogram, MACS2 score/signal summaries, and
  fragment-per-peak counts. These describe the final broad peaks.
- **TSS enrichment:** when `--tss_bed` or `--gtf` is supplied, deepTools builds
  a strand-aware matrix from 3 kb upstream to 3 kb downstream in 10-bp bins
  and emits the matrix and aggregate profile. The current output records the
  computed/skipped status and profile; it does not reduce that curve to a
  single scalar TSS-enrichment score. Without annotation it records
  `skipped_no_annotation`.

These metrics are descriptive QC, not universal pass/fail thresholds. Compare
targets with matched controls and comparable input groups, and establish
project-specific cutoffs before excluding libraries.

## Why broad peaks are primary and narrow peaks are secondary

The primary result reproduces the NanoScope-style, matched-IgG MACS2 BAMPE
call with `--llocal 100000 --keep-dup 1 --broad-cutoff 0.1 --max-gap 1000
--broad`. CUT&Tag signal can occupy extended regulatory domains, and preserving
this broad call is the compatibility objective. Final broad peaks drive
blacklist filtering, FRiP, peak QC, and the principal results.

Narrow calls use the same target/control BAMs, BAMPE mode, local lambda, and
duplicate setting, but exist only to provide localized summits for motif QC.
They never replace the broad results or FRiP. Set
`--motif_use_narrow_peaks false` to skip them and center motif windows on
broad-peak midpoints.

## Motif database and expected motifs

`--motif_db` must be a readable MEME-format DNA motif database. A vertebrate
JASPAR release converted to MEME format, such as
`JASPAR2026_CORE_vertebrates.meme`, is an appropriate starting point provided
its motif names/alternate names match the manifest regular expressions.

For each non-control sample, the pipeline:

1. extracts reference-bounded windows around narrow summits by default, or
   final broad-peak midpoints;
2. creates seeded, non-overlapping, approximately GC-matched genomic
   background that excludes final peaks and an optional blacklist;
3. runs AME Fisher known-motif enrichment, STREME de novo discovery, and FIMO
   scanning;
4. matches the manifest `expected_motif` regular expression and reports motif
   rank, enrichment/effect statistic, significance value, hit count, fraction
   of peaks with a hit, hit-position distribution, and central-hit fraction.

Expected expressions for this dataset are `CTCF`, `GATA1`, and `RUNX1`.
Matching is a regular-expression search against motif ID and alternate ID, so
inspect the chosen database naming before a production run. `motif_not_found`
is an explicit QC failure when a complete AME database report contains no
match; `not_significant` and `no_peaks` are distinct statuses.

For modern AME E-value-only output, the compatibility field named
`best_adjusted_p_value` carries the AME E-value. Interpret it according to the
AME report rather than assuming the field always contains an adjusted
probability.

## Resume behavior

Nextflow stores task state in the work directory and `.nextflow` history. Keep
those paths and rerun the identical command with `-resume`:

```bash
nextflow run main.nf \
  -profile conda \
  -work-dir /shared/nanocut/work \
  --input /shared/project/samples.csv \
  --fasta /shared/references/GRCh38.fa \
  --macs_genome_size hs \
  --outdir /shared/project/results \
  -resume
```

Completed tasks are eligible for cache reuse only when their inputs, command,
configuration, and relevant parameters are unchanged. Moving/changing input
files, deleting the work directory, changing the work path, or changing
parameters can invalidate cache entries. `NXF_CONDA_CACHEDIR` reuses software
environments but is separate from the task cache. Do not delete `work/` until
the run is accepted and no resume is needed.

## Troubleshooting

- **Manifest fails immediately:** run
  `python3 bin/manifest.py validate --input samples.csv --output
  normalized.json` for the precise row-level error. Check the exact header,
  relative FASTQ paths, identifier characters, boolean spelling, repeated
  physical-library paths, and same-input IgG links.
- **FASTQs are unsynchronized:** ensure R1, R2, and I2 contain the same records
  in the same order. `/1` and `/2` suffixes are normalized, but different core
  read IDs or truncated records fail the affected library.
- **Many reads are ambiguous/unassigned:** confirm barcodes are in I2 FASTQ
  orientation and compare observed-I2 counts in demultiplex metrics. Increase
  `--barcode_mismatches` only deliberately; it must remain below barcode
  length, and tied closest barcodes remain ambiguous.
- **Zero-read sample:** correct the manifest/barcode or use
  `--allow_empty true` only for diagnostics. The default failure prevents
  silent downstream empty analyses.
- **Index error:** pass the prefix without `.1.bt2`; verify all six small or all
  six large index files exist and that both sets are not present under the
  same prefix.
- **Motifs do not run:** supply both `--fasta` and a valid `--motif_db`.
  Confirm the MEME header and `MOTIF` records, expected-name regex, usable
  peaks, and background-generation logs.
- **No TSS profile:** supply a valid BED6 `--tss_bed` or a GTF with
  strand-bearing `transcript` features. When both are supplied, BED wins.
- **Conda environment creation is repeated or fails on a cluster:** ensure
  `NXF_CONDA_CACHEDIR` is the same absolute shared path on all nodes, writable
  for creation, and not on node-local storage. Solve each `envs/*.yml`
  independently to diagnose package/channel issues.
- **Docker image failure:** verify Docker access, Linux image architecture,
  registry/network access, and the process-specific image tag. Image pulls
  were not verified here.
- **A resume reruns tasks:** use the same `-work-dir`, preserve `.nextflow`
  history, and compare inputs/parameters/configuration. Review
  `pipeline_info/trace.txt`.

## Testing and limitations

Available checks are:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s tests/unit -p 'test_*.py' -q

for test_script in tests/integration/*.sh; do
  bash -n "$test_script"
  bash "$test_script"
done
```

If `pytest` is installed, the unit suite is also conventionally discoverable
with `pytest tests/unit -q`. Runtime-aware integration scripts execute their
Nextflow harnesses when Nextflow is available and otherwise print an explicit
`SKIP`.

Current limitations:

- Nextflow runtime execution was not verified in this workspace, so no
  successful DSL2, end-to-end, or `-resume` runtime claim is made.
- Conda/Mamba/Micromamba were unavailable, so the pinned environments were
  structurally checked but not solved.
- Real Bowtie2, SAMtools, deepTools, MACS2, BEDTools, MEME Suite, and MultiQC
  execution was not performed here; direct fixtures and offline fakes do not
  replace a production-scale validation.
- Container tags are pinned but were not pulled or architecture-tested.
- The bundled `test` profile uses deterministic fake tools and tiny synthetic
  data; it verifies orchestration contracts, not biological correctness.
- The workflow is specialized to the manifest assay targets `IgG`, `CTCF`,
  `GATA1`, and `RUNX1`.
