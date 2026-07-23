# Bulk nano-CUT&Tag Nextflow Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a manifest-driven Nextflow DSL2 pipeline that demultiplexes synchronized R1/R2/I2 bulk nano-CUT&Tag reads, aligns them, calls IgG-controlled NanoScope-compatible peaks, reports QC, and validates TF motifs.

**Architecture:** Focused Python command-line utilities implement deterministic manifest validation, streaming demultiplexing, peak QC, and motif-summary behavior. Nextflow modules wrap each external bioinformatics tool and connect typed metadata channels; each process uses a pinned, task-specific Conda YAML. Synthetic fixtures support unit and end-to-end test profiles.

**Tech Stack:** Nextflow DSL2, Python 3.12, pytest, Bowtie2, SAMtools, FastQC, deepTools, MACS2, BEDTools, MEME Suite, MultiQC, Conda/Mamba.

## Global Constraints

- Default barcode matching is exact; `--barcode_mismatches` may override it with a non-negative integer smaller than the barcode length.
- Primary peaks use MACS2 BAMPE with `--llocal 100000 --keep-dup 1 --broad-cutoff 0.1 --max-gap 1000 --broad` and the manifest-linked IgG BAM as `-c`.
- Primary FRiP uses final broad peaks and counts paired fragments once.
- Optional narrow peaks are motif-QC inputs only and never replace primary broad peaks or broad-peak FRiP.
- Every process declares a pinned repository-local Conda YAML.
- No production Python behavior is added before its failing test is observed.
- The current workspace does not permit creating `.git`; commit steps are deferred until the project is placed in a writable Git repository.

---

## File Structure

```text
main.nf                         Top-level workflow and channel wiring
nextflow.config                 Profiles, defaults, reports, resources
nextflow_schema.json            Parameter validation and help metadata
modules/local/*.nf              One external-tool operation per module
subworkflows/local/*.nf         Alignment, peak/QC, and motif orchestration
bin/manifest.py                 Manifest validation/normalization CLI
bin/demultiplex_i2.py           Streaming synchronized triple-FASTQ splitter
bin/peak_qc.py                  Fragment FRiP and broadPeak summary CLI
bin/motif_qc.py                 Expected-motif result summarizer
envs/*.yml                      Pinned process environments
assets/samples.example.csv      Complete NX701–NX706 example manifest
conf/test.config                Synthetic end-to-end profile
tests/unit/*.py                 Python unit tests
tests/integration/*.sh          Nextflow smoke/integration checks
tests/data/                     Tiny FASTQs, reference, BED, and motif fixtures
README.md                       Installation, input, execution, outputs, QC
```

### Task 1: Manifest contract and example dataset

**Files:**
- Create: `tests/unit/test_manifest.py`
- Create: `bin/manifest.py`
- Create: `assets/samples.example.csv`

**Interfaces:**
- Produces: `load_and_validate(path: Path) -> list[dict[str, object]]`
- Produces CLI: `manifest.py validate --input samples.csv --output normalized.json`
- Normalized records include `sample_id`, `library_id`, `input_group`, `barcode`, `assay_target`, `is_control`, `control_id`, `expected_motif`, `r1`, `r2`, and `i2`.

- [ ] Write failing pytest cases for the six-library valid manifest, duplicate sample IDs, inconsistent paths within a library, unequal barcode lengths, missing controls, non-IgG controls, and mismatched `input_group`.
- [ ] Run `pytest tests/unit/test_manifest.py -q` and confirm failure because `bin/manifest.py` is absent.
- [ ] Implement CSV parsing, boolean normalization, path resolution relative to the manifest, cross-row validation, and normalized JSON output.
- [ ] Run `pytest tests/unit/test_manifest.py -q` and confirm all cases pass.
- [ ] Add the complete 12-row example manifest mapping 25K/50K/100K targets to NX702/NX703/NX701 IgG respectively, and validate it with the CLI.

### Task 2: Streaming I2 demultiplexer

**Files:**
- Create: `tests/unit/test_demultiplex_i2.py`
- Create: `bin/demultiplex_i2.py`
- Create: `tests/data/demux/` fixtures

**Interfaces:**
- Produces: `hamming(a: str, b: str) -> int`
- Produces: `assign_barcode(observed: str, expected: dict[str, str], max_mismatches: int) -> tuple[str, str | None]`, where status is `assigned`, `ambiguous`, or `unassigned`.
- Produces CLI inputs `--r1`, `--r2`, `--i2`, repeated `--sample SAMPLE_ID=BARCODE`, `--max-mismatches`, `--outdir`, and `--metrics`.

- [ ] Write failing tests for exact A/B assignments, mismatch override, closest-match selection, equal-distance ambiguity, unassigned reads, gzip/plain input, truncated FASTQs, and normalized read-ID mismatches.
- [ ] Run `pytest tests/unit/test_demultiplex_i2.py -q` and observe the missing-module failure.
- [ ] Implement strict four-line FASTQ iteration, gzip detection, identifier normalization, Hamming assignment, and lockstep validation.
- [ ] Run the tests and verify assignment tests pass while output tests fail.
- [ ] Implement gzip R1/R2 writers, atomic metrics JSON/TSV, observed-index counts, empty-sample failure, and `--allow-empty`.
- [ ] Run `pytest tests/unit/test_demultiplex_i2.py -q` and confirm all tests pass.

### Task 3: Peak and FRiP QC utility

**Files:**
- Create: `tests/unit/test_peak_qc.py`
- Create: `bin/peak_qc.py`
- Create: `tests/data/qc/` BED/BEDPE fixtures

**Interfaces:**
- Consumes name-collated BEDPE fragments and broadPeak.
- Produces `summarize_peaks(peaks) -> dict[str, float | int | None]`.
- Produces `calculate_frip(fragments, peaks) -> dict` with `total_fragments`, `fragments_in_peaks`, and `frip`.
- CLI writes per-sample JSON, TSV, width histogram TSV, and fragments-per-peak TSV.

- [ ] Write failing tests proving mates are counted once, a fragment overlapping multiple peaks is counted once for FRiP, empty peaks return FRiP 0, and peak width/quantile/base-coverage summaries are correct.
- [ ] Run `pytest tests/unit/test_peak_qc.py -q` and observe failure because the module is missing.
- [ ] Implement interval validation, merged-overlap FRiP calculation, and peak summaries.
- [ ] Run the focused tests and confirm they pass.
- [ ] Implement the CLI outputs and tests for malformed/empty inputs.
- [ ] Run all unit tests with `pytest tests/unit -q`.

### Task 4: Motif result summarizer

**Files:**
- Create: `tests/unit/test_motif_qc.py`
- Create: `bin/motif_qc.py`
- Create: `tests/data/motifs/` minimal AME/FIMO fixtures

**Interfaces:**
- Produces CLI arguments `--ame`, `--fimo`, `--expected-motif`, `--window`, `--output-prefix`.
- Writes expected-motif JSON/TSV, motif-hit positions TSV, and status values `pass`, `not_significant`, `motif_not_found`, or `no_peaks`.

- [ ] Write failing tests for regex motif selection, adjusted-P ranking, multiple matching motifs, no database match, peak fraction with hits, and centered-position calculations.
- [ ] Run `pytest tests/unit/test_motif_qc.py -q` and observe the expected failure.
- [ ] Implement robust AME/FIMO TSV readers and expected-motif aggregation.
- [ ] Run the focused tests and confirm they pass.
- [ ] Add CLI output tests and run all unit tests.

### Task 5: Conda environments and configuration foundation

**Files:**
- Create: `envs/python.yml`, `envs/fastqc.yml`, `envs/bowtie2.yml`, `envs/samtools.yml`, `envs/deeptools.yml`, `envs/macs2.yml`, `envs/bedtools.yml`, `envs/meme.yml`, `envs/multiqc.yml`
- Create: `nextflow.config`
- Create: `nextflow_schema.json`
- Create: `tests/integration/check_envs.sh`

**Interfaces:**
- Each YAML is accepted by `conda env create --dry-run -f` and contains only direct dependencies.
- Profiles: `conda`, `docker`, and `test`; Conda cache defaults to `${launchDir}/.nextflow-conda` and can be overridden by `NXF_CONDA_CACHEDIR` or user config.

- [ ] Write `check_envs.sh` to assert all expected YAMLs exist, use channel order `conda-forge`, `bioconda`, and contain exact version constraints.
- [ ] Run the checker and observe failure because the YAMLs do not exist.
- [ ] Add minimal pinned YAMLs and rerun the structural checker.
- [ ] Add config defaults for input/reference/output/barcode/MACS/motif options, execution reports, process labels, and profiles.
- [ ] Validate config/schema syntax with installed Nextflow when available; otherwise record the missing-runtime skip distinctly from failure.

### Task 6: Demultiplexing and read-QC Nextflow modules

**Files:**
- Create: `modules/local/validate_manifest.nf`
- Create: `modules/local/demultiplex_i2.nf`
- Create: `modules/local/fastqc.nf`
- Create: `subworkflows/local/demultiplex.nf`
- Create: `tests/integration/test_demultiplex.sh`

**Interfaces:**
- Subworkflow consumes `path manifest` and mismatch parameters.
- Emits derived tuples `tuple(meta, r1, r2)`, metrics files, and versions files.

- [ ] Write a failing synthetic Nextflow test expecting one demultiplex execution per physical library and derived files for both barcodes.
- [ ] Run it and verify failure due to missing modules.
- [ ] Implement manifest validation, parse normalized records into library-grouped channels, and invoke the streaming splitter once per library.
- [ ] Implement derived-sample channel expansion and FastQC.
- [ ] Run the demultiplex integration test and confirm exact assignment metrics and outputs.

### Task 7: Reference, alignment, BAM, and coverage modules

**Files:**
- Create: `modules/local/bowtie2_build.nf`
- Create: `modules/local/bowtie2_align.nf`
- Create: `modules/local/samtools_sort_index.nf`
- Create: `modules/local/samtools_filter.nf`
- Create: `modules/local/samtools_metrics.nf`
- Create: `modules/local/bamcoverage.nf`
- Create: `subworkflows/local/align_qc.nf`
- Create: `tests/integration/test_alignment.sh`

**Interfaces:**
- Consumes derived tuples plus FASTA/index.
- Emits metadata-keyed analysis BAM, filtered BAM, indexes, bigWig, alignment metrics, and versions.

- [ ] Write a failing integration test using the miniature reference and deterministic paired reads.
- [ ] Implement index construction and Bowtie2 paired alignment with captured summary.
- [ ] Implement coordinate sorting/indexing and filtered properly paired primary BAM at `--min_mapq` default 5.
- [ ] Add flagstat/stats/idxstats and duplicate/library-complexity metrics.
- [ ] Add NanoScope-style `bamCoverage` settings.
- [ ] Run alignment integration tests and assert mapped pairs, indexed BAMs, metrics, and bigWig exist.

### Task 8: IgG-controlled broad and optional narrow peaks

**Files:**
- Create: `modules/local/macs2_broad.nf`
- Create: `modules/local/macs2_narrow.nf`
- Create: `modules/local/filter_blacklist.nf`
- Create: `subworkflows/local/peaks.nf`
- Create: `tests/integration/test_peaks.sh`

**Interfaces:**
- Joins target metadata to control BAM by `control_id`.
- Emits primary raw/final broadPeak and optional secondary narrowPeak/summits.

- [ ] Write a failing test that inspects task commands and requires `-c`, `-f BAMPE`, `--llocal 100000`, `--keep-dup 1`, `--broad-cutoff 0.1`, `--max-gap 1000`, and `--broad`.
- [ ] Implement target/control keyed joining with explicit failure for unmatched controls.
- [ ] Implement primary broad calls and confirm IgG is not self-called.
- [ ] Implement optional blacklist filtering while retaining raw peaks.
- [ ] Implement optional controlled narrow calls used only by motif channels.
- [ ] Run peak integration tests and verify target/control pairing for all input groups.

### Task 9: Pipeline QC processes and report aggregation

**Files:**
- Create: `modules/local/bam_to_fragments.nf`
- Create: `modules/local/peak_qc.nf`
- Create: `modules/local/tss_enrichment.nf`
- Create: `modules/local/multiqc.nf`
- Create: `subworkflows/local/qc.nf`
- Create: `tests/integration/test_qc.sh`

**Interfaces:**
- Emits target QC JSON/TSV/plots, combined summary TSV, optional TSS profiles, and MultiQC HTML/data.

- [ ] Write a failing test requiring fragment-based FRiP, peak-width outputs, combined target summary, and graceful empty-peak handling.
- [ ] Implement BAM-to-fragment conversion with one record per proper pair.
- [ ] Wrap `peak_qc.py` and publish per-target outputs.
- [ ] Implement optional strand-aware TSS BED creation and deepTools enrichment.
- [ ] Create MultiQC custom-content inputs for demultiplexing, FRiP, peak count, expected motif, and library QC.
- [ ] Run QC integration tests and verify IgG has library QC but no default target FRiP.

### Task 10: Motif workflow

**Files:**
- Create: `modules/local/prepare_motif_sequences.nf`
- Create: `modules/local/ame.nf`
- Create: `modules/local/streme.nf`
- Create: `modules/local/fimo.nf`
- Create: `modules/local/motif_summary.nf`
- Create: `subworkflows/local/motifs.nf`
- Create: `tests/integration/test_motifs.sh`

**Interfaces:**
- Consumes final motif peaks/summits, FASTA, optional blacklist, motif DB, and `expected_motif` metadata.
- Emits known/de novo motif directories, scans, and expected-motif QC summaries.

- [ ] Write a failing fixture test requiring bounded 200-bp windows, shuffled non-overlapping background, AME/STREME/FIMO outputs, and explicit motif-not-found status.
- [ ] Implement summit-or-midpoint selection and reference-bound clipping.
- [ ] Implement reproducible shuffled background with fixed seed, peak exclusion, and optional blacklist exclusion.
- [ ] Wrap AME, STREME, and FIMO with the MEME environment.
- [ ] Wrap `motif_qc.py`, publish expected-motif metrics, and feed them to MultiQC custom content.
- [ ] Run motif integration tests with the miniature MEME database.

### Task 11: Top-level workflow and synthetic end-to-end profile

**Files:**
- Create: `main.nf`
- Create: `conf/test.config`
- Create: `tests/data/e2e/` fixtures
- Create: `tests/integration/test_e2e.sh`

**Interfaces:**
- `nextflow run main.nf -profile test` is offline, deterministic, and exercises demultiplexing through reporting.
- Production invocation accepts all schema-documented parameters.

- [ ] Write the failing end-to-end test with expected output paths and summary values.
- [ ] Wire subworkflows in dependency order and merge all versions/metrics channels.
- [ ] Add workflow-level parameter/path validation and completion summary.
- [ ] Add synthetic reference, FASTQs/I2, blacklist, TSS, and motif database fixtures.
- [ ] Run `nextflow run main.nf -profile test` and fix only implementation defects until it succeeds.
- [ ] Rerun with `-resume` and verify completed processes are cached.

### Task 12: Documentation and final verification

**Files:**
- Create: `README.md`
- Create: `.gitignore`
- Create: `CHANGELOG.md`
- Create: `outputs/` packaged pipeline archive only after verification

**Interfaces:**
- README covers prerequisites, Conda cache configuration, manifest semantics, six-library mapping, reference inputs, execution, outputs, QC interpretation, motif interpretation, resume behavior, and troubleshooting.

- [ ] Write documentation assertions checking every CLI parameter and environment file is mentioned or discoverable in generated help.
- [ ] Write README examples for local Conda and shared-HPC cache execution.
- [ ] Run `pytest tests/unit -q`.
- [ ] Run every `tests/integration/test_*.sh` test.
- [ ] Run `nextflow run main.nf -profile test` from a clean work directory.
- [ ] Run the verification-before-completion checklist, inspect generated reports, and record exact commands/results.
- [ ] Package the verified pipeline into `outputs/bulk-nanocut-tag-nextflow.tar.gz` without work/cache directories.

