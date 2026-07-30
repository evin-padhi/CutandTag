# Consolidated QC Dashboard Design

## Purpose

The bulk nano-CUT&Tag workflow currently produces useful per-process QC files and
a MultiQC report, but users must move between several directories to compare
samples. This change will add a consolidated, self-contained report and
machine-readable run summaries while preserving the existing MultiQC output.

The report is descriptive. It will expose missing or failed QC inputs, but it
will not assign biological pass/fail labels or impose fixed QC thresholds.

## User-facing outputs

The workflow will continue to publish:

```text
<outdir>/reports/multiqc/multiqc_report.html
```

It will additionally publish:

```text
<outdir>/reports/qc_dashboard/qc_dashboard.html
<outdir>/reports/qc_dashboard/qc_summary.tsv
<outdir>/reports/qc_dashboard/qc_summary.json
<outdir>/reports/qc_dashboard/top_motifs.tsv
<outdir>/reports/qc_dashboard/tss_profiles.tsv
```

`qc_dashboard.html` will be self-contained so it can be copied from the HPC and
opened locally without a web server or network connection.

The existing per-sample FRiP source remains:

```text
<outdir>/qc/peaks/<sample_id>/<sample_id>.peak_qc.tsv
<outdir>/qc/peaks/<sample_id>/<sample_id>.peak_qc.json
```

The TSV stores FRiP in the row named `frip`.

## Report organization

### Run overview

The first section will show:

- sample count and target/control composition;
- sample, antibody/target, and control metadata;
- availability of each expected QC family;
- explicit warnings for missing, empty, skipped, or failed optional QC outputs.

Missing optional values will be rendered as `NA`. A missing optional metric must
not make report generation fail. Structurally invalid inputs that could silently
mislabel samples will remain fatal and produce an actionable error.

### Demultiplexing

For each parent library, the report will show:

- total input reads;
- assigned reads;
- unassigned reads;
- assigned percentage;
- per-derived-sample read allocation where available.

These values will also be included in the joined per-sample summary.

### Alignment and library QC

All target libraries and IgG controls will be shown together for technical QC:

- raw and retained read counts;
- mapped percentage;
- properly paired percentage;
- MAPQ-filtered reads, fragments, and retained fraction;
- duplicate percentage;
- mitochondrial percentage;
- estimated library size;
- insert-size summary and distribution where available.

Plots will retain a visible distinction between target libraries and controls.

### Peak calling and FRiP

The peak section will show:

- peak count;
- total bases covered by peaks;
- peak-width summary and distribution;
- usable fragment count;
- fragments overlapping peaks;
- FRiP.

IgG controls will remain visible but will be visually separated from target
libraries so that controls are not interpreted as equivalent biological
enrichment samples. The report will not color samples as passing or failing.

### TSS enrichment

The report will include:

- overlaid TSS profiles for cross-sample comparison;
- a per-sample scalar TSS enrichment score;
- the existing computation/skipped status.

The scalar score will use the normalized aggregate profile and follow a
transparent center-to-flank definition:

```text
TSS enrichment = center signal / mean flank signal
```

The center signal is the signal at the profile center. The flank signal is the
mean of the terminal 100 bp at both ends of the profile. If the profile has no
finite flank signal, or the flank mean is zero, the score will be `NA` and a
warning will be recorded. The report generator will record the exact window
definition in the report and JSON metadata.

This complements, rather than replaces, visual inspection of the full profile.

### Motif enrichment

For each target sample, the report will show:

- expected-motif QC status;
- best matching expected motif identifier;
- expected-motif adjusted significance statistic;
- the top ten AME motifs ranked by AME significance;
- the AME computed/skipped/failed status.

IgG samples will be identified as controls and excluded from biological
expected-motif interpretation. If motif outputs are absent because the sample
has no peaks or the task was skipped, the report will state that explicitly.

## Data model

`qc_summary.tsv` will contain one row per derived sample. It will be the main
spreadsheet/R-friendly table and will join sample metadata with demultiplexing,
alignment, library, peak, FRiP, TSS, and expected-motif metrics.

`qc_summary.json` will contain:

- schema and report-generator versions;
- run-level counts and warnings;
- metric definitions, including the TSS calculation;
- one structured object per sample;
- input availability/status information.

`top_motifs.tsv` will contain up to ten AME rows per target sample, with rank,
motif identifiers, significance values, and available enrichment counts.

`tss_profiles.tsv` will use a tidy layout with one row per sample and relative
TSS position so it can be replotted independently.

Stable, documented column names will be used. Numeric missing values will be
empty in TSV and `null` in JSON, not ambiguous strings or zeroes.

## Workflow architecture

The existing MultiQC process remains the quick overview and compatibility
report. A new local report process will consume the already generated custom QC
tables plus the detailed peak, TSS, and motif artifacts needed for plots.

A repository-owned Python report generator will:

1. validate and parse the input QC artifacts;
2. join records by sample identifier;
3. calculate the TSS scalar;
4. write the four machine-readable/HTML outputs;
5. embed report styling, plot data, and plotting code into one HTML file.

The report generator will use the workflow's managed environment. It will avoid
external CDNs and runtime downloads. The Nextflow process will publish the
report directory and expose the HTML and joined summary as workflow outputs.

MultiQC will gain a clear link or pointer to the detailed report where the
supported MultiQC custom-content mechanism permits it. The workflow completion
summary and README will list both report locations.

## Error handling

- Missing optional per-sample metrics become `NA` plus visible warnings.
- Empty or skipped peak/motif/TSS results retain their explicit status.
- Duplicate sample identifiers, incompatible schemas, or contradictory sample
  metadata are fatal.
- Non-finite numeric values are rejected or converted to missing with a warning,
  depending on whether they can corrupt joins.
- Report generation writes outputs atomically so a failed process does not leave
  a plausible but incomplete HTML report.

## Verification

Unit tests will cover:

- parsing each supported QC input;
- joined per-sample output and stable columns;
- demultiplexing percentages;
- center-to-flank TSS enrichment;
- zero/missing flank behavior;
- AME top-ten ranking;
- target/control classification;
- incomplete and skipped samples;
- deterministic JSON, TSV, and HTML generation.

Integration tests will cover:

- Nextflow channel wiring;
- publication of all five new report outputs;
- a mixed target/IgG fixture;
- successful report generation when optional QC data are missing;
- preservation of the existing MultiQC report and combined peak summary.

Documentation will describe the report paths, metric definitions, and the
location of the original per-sample FRiP files.

## Non-goals for this pull request

- automatic biological pass/fail thresholds;
- cohort-specific acceptance criteria;
- statistical comparison between experimental conditions;
- differential binding or differential peak analysis;
- replacing MultiQC;
- requiring an interactive server after workflow completion.
