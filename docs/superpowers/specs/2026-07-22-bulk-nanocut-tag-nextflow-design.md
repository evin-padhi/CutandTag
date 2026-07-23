# Bulk nano-CUT&Tag Nextflow Pipeline Design

## Objective

Build a reproducible Nextflow DSL2 pipeline that demultiplexes bulk nano-CUT&Tag paired-end reads using an 8-base I2 modality barcode, aligns each derived sample to a reference genome, calls NanoScope-compatible broad peaks with matched IgG controls, calculates library and peak QC, and validates transcription-factor experiments by motif enrichment.

The initial dataset contains six NextSeq 500 libraries and approximately 60 million read pairs in total:

| Input | IgG + CTCF library | GATA1 + RUNX1 library |
|---|---|---|
| 25K PBMCs | NX702 | NX704 |
| 50K PBMCs | NX703 | NX705 |
| 100K PBMCs | NX701 | NX706 |

Each read pair has one synchronized I2 record. Barcode A is `TATAGCCT`; barcode B is `ATAGAGGC`. Within NX701–NX703, barcode A identifies IgG and barcode B identifies CTCF. Within NX704–NX706, barcode A identifies GATA1 and barcode B identifies RUNX1.

## User Interface

The pipeline is launched with a manifest, reference FASTA, output directory, and optional annotations:

```bash
nextflow run main.nf \
  --input samples.csv \
  --fasta GRCh38.fa \
  --outdir results \
  --blacklist hg38-blacklist.v2.bed \
  --gtf gencode.annotation.gtf \
  --motif_db JASPAR2026_CORE_vertebrates.meme \
  -profile docker
```

Conda and Docker profiles will be supplied. Tool versions will be pinned. A local profile will permit execution without a scheduler, and resource settings will be overridable through a user configuration file.

### Process environments

Every Nextflow process will declare a repository-local Conda environment with the `conda` directive. Environments are split by responsibility so a task installs only the tools it uses and dependency conflicts remain isolated:

```text
envs/
  python.yml       # manifest validation, demultiplexing, and custom QC utilities
  fastqc.yml       # FastQC
  bowtie2.yml      # Bowtie2 and reference-index construction
  samtools.yml     # BAM conversion, sorting, filtering, indexing, and metrics
  deeptools.yml    # bamCoverage and optional TSS enrichment
  macs2.yml        # MACS2 broad and narrow peak calling
  bedtools.yml     # blacklist filtering, interval operations, and sequence windows
  meme.yml         # AME, STREME, FIMO, and motif utilities
  multiqc.yml      # MultiQC report assembly
```

Processes that perform the same operation share the corresponding environment; unrelated tools are not combined into a monolithic environment. Each YAML will use strict channel order (`conda-forge`, then `bioconda`, then `defaults` only when required), pin exact package versions, and include only direct runtime dependencies. The `conda` profile enables `conda.enabled = true`; users may select Mamba/Micromamba through Nextflow configuration for faster environment creation. Nextflow's environment cache remains configurable through `conda.cacheDir`, allowing cluster nodes to reuse environments across tasks and resumed runs.

The Docker profile will remain available as an alternative and will use process-specific BioContainers-compatible images wherever practical. Conda YAML files are the canonical dependency definitions and are exercised by automated environment-creation checks.

## Manifest

The CSV contains one row per expected barcode-derived sample. Repeating the three FASTQ paths is intentional: the pipeline groups rows by `library_id` and demultiplexes each physical library exactly once.

Required columns:

- `sample_id`: unique identifier for the derived sample.
- `library_id`: identifier shared by rows originating from the same physical library.
- `input_group`: biological input or grouping label, such as `25K`.
- `barcode`: expected I2 sequence in FASTQ orientation.
- `assay_target`: `IgG`, `CTCF`, `GATA1`, or `RUNX1` for this dataset.
- `is_control`: boolean indicating an IgG control.
- `control_id`: `sample_id` of the matched IgG; blank only for controls.
- `expected_motif`: motif-name regular expression used to locate the expected motif in database results; blank for controls.
- `r1`, `r2`, `i2`: synchronized gzip-compressed or uncompressed FASTQ paths.

Example rows:

```csv
sample_id,library_id,input_group,barcode,assay_target,is_control,control_id,expected_motif,r1,r2,i2
NX701_IgG,NX701,100K,TATAGCCT,IgG,true,,,data/NX701_R1.fastq.gz,data/NX701_R2.fastq.gz,data/NX701_I2.fastq.gz
NX701_CTCF,NX701,100K,ATAGAGGC,CTCF,false,NX701_IgG,CTCF,data/NX701_R1.fastq.gz,data/NX701_R2.fastq.gz,data/NX701_I2.fastq.gz
NX706_GATA1,NX706,100K,TATAGCCT,GATA1,false,NX701_IgG,GATA1,data/NX706_R1.fastq.gz,data/NX706_R2.fastq.gz,data/NX706_I2.fastq.gz
NX706_RUNX1,NX706,100K,ATAGAGGC,RUNX1,false,NX701_IgG,RUNX1,data/NX706_R1.fastq.gz,data/NX706_R2.fastq.gz,data/NX706_I2.fastq.gz
```

Validation fails before computation when sample IDs are duplicated, barcodes are not equal length, a library points to inconsistent FASTQ paths, a target references a missing/non-control control, a control references another control, files are missing, or required values are blank. A target may use an IgG from a different library, but their `input_group` values must match.

## Data Flow

### 1. Manifest validation and grouping

A small Python utility validates the complete manifest and emits normalized JSON records. Nextflow creates one demultiplexing task per `library_id` and one downstream task per `sample_id`.

### 2. I2 demultiplexing

A streaming Python utility reads R1, R2, and I2 in lockstep without loading the library into memory. It verifies record counts and normalized read identifiers. For each I2 sequence, it computes Hamming distance to every expected barcode for that library.

- `--barcode_mismatches` defaults to `0` and must be a non-negative integer smaller than the barcode length.
- A unique closest barcode within the threshold is assigned.
- A tie within the threshold is classified as ambiguous.
- No barcode within the threshold is classified as unassigned.
- Assigned R1 and R2 records are written to gzip-compressed per-sample FASTQs.
- I2 is not propagated downstream.

The task emits per-barcode read counts, assigned/ambiguous/unassigned counts and fractions, and observed-I2 sequence counts. Empty derived samples are an error by default, with an explicit option to permit them for diagnostic runs.

### 3. Read QC and alignment

FastQC runs on demultiplexed R1/R2. Reads are aligned with Bowtie2 in paired-end mode against an index built from `--fasta` or a supplied `--bowtie2_index`. SAMtools converts, coordinate-sorts, indexes, and collects `flagstat`, `stats`, insert-size, and idxstats-derived mitochondrial metrics.

The primary analysis BAM retains all reads so that MACS2 can reproduce NanoScope duplicate handling. A filtered BAM containing properly paired primary alignments above configurable mapping quality (default MAPQ 5) is used for coverage and QC. Duplicate metrics are reported; removal is not applied before primary peak calling because NanoScope passes `--keep-dup 1` to MACS2.

### 4. Coverage

deepTools `bamCoverage` reproduces NanoScope settings on the filtered BAM:

```text
--minMappingQuality 5 --binSize 50 --centerReads --smoothLength 250
--normalizeUsing RPKM --ignoreDuplicates --extendReads
```

Each sample receives an indexed BAM, BAM index, and RPKM bigWig.

### 5. Primary peak calling

Every non-control sample is called with MACS2 against the manifest-linked IgG BAM. The primary broad-peak command preserves NanoScope parameters and adds the control:

```text
macs2 callpeak -t TARGET.bam -c IGG.bam -g GENOME_SIZE -f BAMPE \
  -n SAMPLE --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 \
  --max-gap 1000 --broad
```

`--macs_genome_size` is required and accepts a MACS shortcut such as `hs` or an effective genome-size integer. IgG samples are not peak-called against themselves. MACS2 logs and all peak outputs are retained.

If a blacklist is supplied, `bedtools intersect -v` produces a filtered broadPeak file. Both original and filtered peaks remain available. Filtered peaks are used for final peak QC, FRiP, and motif analysis.

### 6. Optional narrow peaks for motif QC

When `--motif_use_narrow_peaks true` (the default), each non-control sample also receives an IgG-controlled MACS2 narrow call using BAMPE, `--llocal 100000`, and `--keep-dup 1`. These narrow peaks are explicitly secondary and are used only for motif analysis. They do not replace broad peaks in the principal results or FRiP calculation. Setting the parameter to false makes motif analysis use broad-peak midpoints.

### 7. Library and peak QC

The pipeline reports, per sample:

- raw and assigned read-pair counts;
- barcode assignment, ambiguity, and unassigned fractions;
- FastQC results;
- overall alignment rate, properly paired fraction, MAPQ-filtered fraction, duplicate fraction, mitochondrial fraction, insert-size distribution, and library-complexity metrics;
- peak count, total bases covered, genome coverage, peak-width minimum/mean/median/maximum and quantiles;
- peak-width histogram and cumulative distribution;
- MACS2 score, signal-value, and reads/fragments-per-peak distributions;
- number and fraction of original peaks removed by an optional blacklist;
- FRiP.

FRiP is defined as the number of unique properly paired fragments in the filtered BAM overlapping at least one final broad peak, divided by all unique properly paired fragments in that BAM. A fragment is counted once even when both mates or multiple peaks overlap. IgG libraries receive library QC but do not receive target FRiP by default.

When `--gtf` or `--tss_bed` is supplied, deepTools calculates TSS enrichment. The GTF path is converted to strand-aware TSS positions. TSS enrichment is skipped, with a clear report entry, when neither annotation is supplied.

Per-sample TSV/JSON metrics, combined summary tables, plots, and MultiQC inputs are retained. MultiQC provides the top-level run report.

### 8. Motif validation

Motif analysis runs for non-control samples when `--motif_db` points to a MEME-format database.

- Narrow-peak summits are used by default. When narrow analysis is disabled, broad-peak midpoints are used.
- Fixed-width, reference-bounded windows are extracted around summits/midpoints; `--motif_window` defaults to 200 bp total width.
- MEME Suite AME tests known-motif enrichment against length- and approximately GC-matched shuffled genomic background that excludes called peaks and an optional blacklist.
- STREME performs de novo motif discovery against the same background.
- FIMO scans peak windows for database motifs.
- The manifest `expected_motif` regular expression selects one or more matching database motifs. No match is a reported QC failure for that sample, not a silent skip.
- Reports include expected-motif enrichment rank, effect/enrichment statistic, adjusted P value, fraction of peaks with the expected motif, hit counts, motif-position distribution, and a central-enrichment summary.

IgG sequences are not the sole statistical background because IgG can yield too few reliable regions. Where an IgG peak set is available, the pipeline may report it as an additional descriptive comparison while retaining matched genomic background for the primary test.

## Output Layout

```text
results/
  demultiplex/<library_id>/
  fastqc/<sample_id>/
  alignment/<sample_id>/
  coverage/<sample_id>/
  peaks/<sample_id>/broad/
  peaks/<sample_id>/narrow_motif_qc/
  qc/library/
  qc/peaks/
  motifs/<sample_id>/ame/
  motifs/<sample_id>/streme/
  motifs/<sample_id>/fimo/
  reports/multiqc/
  reports/summary/
  pipeline_info/
```

The `pipeline_info` directory includes validated parameters, software versions, execution reports, timeline, trace, DAG, and the normalized manifest.

## Failure Behavior

- Manifest problems stop the run before expensive processes begin.
- Unsynchronized FASTQs stop only the affected library task with the offending record number and identifiers.
- Ambiguous and unassigned reads are counted but not aligned.
- A zero-read derived sample fails by default.
- A target with no usable matched IgG BAM cannot be peak-called and causes a clear failure.
- Empty peak sets remain valid outputs and generate zero/NA QC values; motif analysis is skipped with a recorded reason.
- Missing optional blacklist, annotation, or motif database features are skipped only when the corresponding parameter is absent, never when an invalid path was supplied.

## Testing Strategy

Unit tests use `pytest` and tiny synthetic FASTQs/reference data. Tests are written before implementation and cover:

- manifest validation and matched-control relationships;
- exact demultiplexing;
- configurable mismatch assignment;
- tie-to-ambiguous behavior;
- unassigned reads;
- gzip and plain-text FASTQs;
- truncated files and mismatched read identifiers;
- fragment-based FRiP without double-counting mates or overlapping peaks;
- peak-width and peak-summary calculations;
- expected-motif lookup and summary parsing.

Nextflow tests use nf-test or deterministic smoke-test scripts to verify channel construction, one demultiplex task per physical library, correct target/IgG pairing, MACS2 broad command parameters, optional narrow motif calls, and expected publication paths. A miniature end-to-end profile runs entirely on bundled synthetic data and emits non-empty alignment, broad-peak, QC, and report artifacts. Motif tools are covered by a small dedicated integration fixture rather than requiring the full JASPAR database.

## Acceptance Criteria

The pipeline is complete when:

1. The supplied six-library manifest validates and maps each target to the same-input IgG.
2. Every physical library is streamed once and split correctly by `TATAGCCT` and `ATAGAGGC` with exact matching by default.
3. A resumed run does not repeat completed processes when inputs and parameters are unchanged.
4. Every derived sample has demultiplexing, FastQC, alignment, BAM, index, coverage, and library-QC outputs.
5. Every target has an IgG-controlled NanoScope-compatible broadPeak file and FRiP/peak-QC outputs.
6. Every TF sample has known and de novo motif results when a motif database is supplied, including an explicit expected-motif QC result.
7. The synthetic end-to-end run and all unit tests pass under a documented supported profile.
8. The README documents installation, manifest construction, reference preparation, execution, outputs, interpretation, and troubleshooting.
