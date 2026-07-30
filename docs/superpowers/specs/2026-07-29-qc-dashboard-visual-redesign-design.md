# QC Dashboard Visual Redesign

## Purpose

Redesign the consolidated QC HTML report so that its plots can be interpreted
without repeatedly consulting the raw tables. The report will follow the
publication-style visual language established by the supplied sequencing,
peak, TSS, and motif figures while remaining a deterministic, self-contained
HTML file that can be copied from an HPC and opened offline.

The redesign will:

- label every quantitative axis with a metric and unit;
- put compact values on bars;
- use consistent antibody-target colors across the report;
- compare samples with normalized distributions rather than raw counts;
- replace per-base-pair HTML tables with compact plotted summaries;
- add fragments-per-peak coverage distributions;
- add a target-aware motif-enrichment heatmap;
- retain full-resolution machine-readable data outside the visible HTML.

The report remains descriptive. It will not classify a sample as passing or
failing or introduce biological thresholds.

## Selected Approach

Use a publication-style static dashboard rendered with inline SVG and CSS.
Small hover titles may supplement visible labels, but no plot will depend on
JavaScript, an external font, a network request, or an external plotting
library. This keeps the report reproducible, printable, accessible, and easy
to run in the current Python-only dashboard environment.

## Report Structure

The report will contain these sections:

1. Run overview and input-family availability.
2. Sequencing, demultiplexing, and alignment QC.
3. Peak calling and peak-level QC.
4. Insert-size, peak-width, and fragments-per-peak distributions.
5. TSS profiles and enrichment scores.
6. Motif-enrichment heatmap and expected-motif summary.
7. Warnings.
8. Collapsed supporting data tables.

The run overview and warnings remain visible. Detailed tables use accessible
`details`/`summary` controls and are collapsed by default.

## Visual Language

### Target colors

The report will assign colors by assay target rather than by sample:

- IgG control: neutral grey;
- CTCF: blue;
- GATA1: green;
- RUNX1: orange.

Additional targets receive deterministic colors from an accessible palette.
The same target always receives the same color in every panel. Controls retain
a secondary distinction through a hatch, square marker, or dashed stroke so
the report does not depend on color alone.

Individual samples within one target use deterministic line styles or markers.
Overlaid line plots end in direct sample labels. A static label-packing
algorithm will sort endpoints, enforce a minimum vertical separation, constrain
labels to the plot area, and draw leader lines back to displaced endpoints.
A compact legend remains available as an accessibility and dense-cohort
fallback.

### Units and value labels

Every quantitative plot will show an explicit axis title and readable ticks.
Bars will carry visible values above or beside them. Formatting depends on the
metric:

- read and fragment counts: millions (`M`) where appropriate;
- peak counts: thousands (`K`);
- covered peak bases: megabases (`Mb`);
- percentages and FRiP: percent;
- widths: base pairs or kilobases;
- significance: `-log10(adjusted p)`.

Tooltips preserve the unshortened value. A transformed axis must say so in its
title; the report will never silently plot a logarithm.

Numeric values in visible HTML tables will use three significant digits.
Large counts may use the same `K`, `M`, or `Mb` suffixes as their plots, and
very small values may use scientific notation. Sample IDs, motif IDs, ranks,
integer bin boundaries, and categorical fields remain exact. This is a
presentation-only transformation: TSV and JSON outputs retain their existing
full numeric precision.

## Sequencing and Alignment Panels

The sequencing section will use compact small multiples for:

- assigned read pairs in millions;
- within-library barcode balance as percent of assigned reads;
- mapped reads as percent;
- usable MAPQ-filtered fragments in millions on a log10 axis;
- PCR duplicate rate as percent;
- end-to-end usable yield as percent of assigned reads.

The charts will share target colors and sample ordering. The metric definition
and denominator will appear beside each panel or in a concise section note.
Undefined ratios will be shown as `NA`, not zero.

## Peak Panels

The peak section will show:

- peak count in thousands;
- FRiP as percent;
- total bases covered by peaks in megabases;
- peak-width median with Q25–Q75 and minimum–maximum ranges;
- peak count versus usable fragments, with both axes labeled and samples
  directly labeled.

IgG controls have no pipeline peak calls and will remain explicitly
not-applicable rather than appearing as zero-valued target samples.

## Distribution Transformations

All distribution plots use a shared target palette and one line per sample.
They compare distribution shape rather than raw library depth.

### Insert size

Insert sizes will be aggregated into shared, left-closed 250-bp bins:

```text
[0, 250), [250, 500), [500, 750), ...
```

Each bin is normalized to percent of the sample's observed insert-size pairs.
All samples share every bin from zero through the cohort-wide observed maximum,
including zero-valued bins. The x-axis is linear and labeled in base pairs.
The binned HTML table contains sample ID, bin start, bin end, pair count, and
within-sample percent.

### Peak width

Peak widths will be aggregated into the same style of shared, left-closed
250-bp bins:

```text
[0, 250), [250, 500), [500, 750), ...
```

Each bin is normalized to percent of the sample's called peaks. All target
samples share every bin from zero through the cohort-wide observed maximum,
including zero-valued bins. The x-axis is linear and labels values in bp/kb.
The binned HTML table contains sample ID, bin start, bin end, peak count, and
within-sample percent.

### Fragments per peak

The existing `*.peak_qc.fragments_per_peak.tsv` output is the raw source.
The dashboard will render an ECDF of fragment count per called peak with a
logarithmic x-axis and cumulative-percent y-axis. Zero-fragment peaks are
counted, reported explicitly, and drawn at a labeled zero origin before the
positive logarithmic domain.

The dashboard process will reduce the raw per-peak input to compact ECDF data;
it will not embed every peak record in the HTML.

## TSS Panels

TSS presentation will retain the aggregate profiles and scalar enrichment
scores. Profiles use the target palette, sample-specific strokes, direct
endpoint labels, and a visible TSS reference line. Axes will read:

- x: distance from TSS (bp);
- y: mean coverage (RPKM).

The scalar TSS enrichment chart will be sorted by score, colored by target,
and show the numeric score beside each bar. The existing center-to-terminal
flank definition remains unchanged and visible in the report.

## Motif Heatmap

### Matrix

Heatmap columns are target samples. Rows include:

1. every motif cognate to an expected assay TF;
2. the cohort-wide top 15 remaining motifs, ranked by strongest adjusted
   significance observed in any sample.

Duplicate motif identities are collapsed deterministically using normalized
motif ID plus alternate name. Cognate rows appear first, grouped by expected
TF, followed by remaining rows ordered by significance and stable lexical
tie-breakers.

Each cell represents `-log10(adjusted p)`, capped at 60. The cap is displayed
in the legend and cells at the cap use a `>60` label. Zero adjusted p-values
are assigned to the cap. Adjusted p-values greater than `0.1`, missing motif
results, and unavailable comparisons display `ns` with distinct accessible
titles that preserve the underlying reason.

### Cognate matching

Cognate matching is case-insensitive and token-aware across motif ID and
alternate name. Tokens are delimited by non-alphanumeric characters. Expected
`GATA1` therefore matches:

- `GATA1`;
- `GATA1::TAL1`;
- `TAL1::GATA1`.

It does not match `GATA10`. The same rule applies to CTCF, RUNX1, and future
targets. Every cognate cell receives a visible outline, including complex
motifs rather than only exact-name matches.

The collapsed supporting table contains the sample, motif ID, alternate name,
rank, adjusted p-value, transformed significance, and cognate flag for every
displayed heatmap cell.

## Data Flow and Contracts

The current full-resolution insert-size and peak-width arrays remain unchanged
in `qc_summary.json`. Their visible HTML tables are replaced by compact plotted
summaries.

Heatmap construction will consume the complete AME tables before the existing
top-ten export truncation. `top_motifs.tsv` retains its current bounded
consumer contract.

The peak-QC process already emits per-peak fragment counts. The QC subworkflow
will carry those files into `QC_DASHBOARD` as a separately staged, sample-keyed
input. The dashboard parser will validate:

- exactly one staged file per expected target sample;
- unique sample identities;
- non-negative integer fragment counts;
- valid headers;
- agreement between filename, metadata, and staged sample order.

The public JSON gains an additive
`peak.fragments_per_peak_distribution` histogram whose rows contain
`fragment_count` and `peak_count`. The dashboard generator version increments
from `1.0.0` to `1.1.0` to document the additive report shape. Existing scalar
columns and standalone raw peak-QC files remain unchanged.

## Component Boundaries

The implementation will keep data transforms separate from SVG rendering:

- target color and sample-style assignment;
- compact number and tick formatting;
- three-significant-digit HTML table formatting;
- bar-chart geometry and value labeling;
- shared 250-bp insert-size and peak-width binning and normalization;
- ECDF construction and downsampling;
- endpoint-label packing;
- motif identity normalization and cognate matching;
- cohort motif selection and heatmap matrix construction;
- reusable SVG axes, bars, lines, markers, ranges, scatter points, and cells;
- collapsed table rendering.

Each transform accepts validated Python data and returns a deterministic,
plain-data representation. SVG helpers consume those representations without
recomputing biological metrics.

## Missing and Malformed Data

Missing optional input produces an unavailable panel with a concise reason and
does not abort the workflow. Empty observed inputs remain distinct from missing
inputs. Controls remain not-applicable for peak and motif analyses.

Present but malformed input remains fatal. Errors will identify the file and
the violated contract, including invalid numeric values, duplicate samples,
noncontiguous staged ordinals, or mismatched metadata. A partial dashboard
bundle will not be published.

## Accessibility and Scaling

Every SVG will have a descriptive accessible name. Tables retain unique region
labels. Color is supplemented by labels, markers, hatches, or strokes. Text
and leader lines are part of the SVG so printed and offline copies retain
sample identities.

Bar panels may scroll horizontally for large cohorts. Distribution and TSS
plots remain fixed-width and use label packing plus a legend. Heatmap cells
have a minimum readable size and the heatmap scrolls within its section rather
than expanding the page width.

## Verification

Unit tests will cover:

- exact 250-bp insert-size and peak-width bin boundaries and within-sample
  percentages;
- ECDF monotonicity, endpoints, repeated values, empty inputs, and zero counts;
- compact-number and tick formatting;
- three-significant-digit HTML table values without loss of TSV/JSON
  precision;
- target-color determinism and control distinctions;
- endpoint-label separation and plot-bound constraints;
- token-aware cognate matches and false-positive exclusions;
- top-15 plus forced-cognate motif selection;
- heatmap significance clipping and `ns` states;
- visible axis titles, units, and bar-value labels;
- collapsed compact tables and preservation of full-resolution JSON.

Integration tests will cover:

- Nextflow transport and staging of fragments-per-peak files;
- sample identity and count validation at the dashboard boundary;
- complete output publication and resume-safe wiring;
- a large synthetic cohort without runaway HTML size or overlapping bar
  labels.

A representative dashboard will be rendered and visually inspected at desktop
and narrow viewport widths. Visual QA will check label collisions, axis
legibility, heatmap scrolling, table collapse behavior, target-color
consistency, and print-friendly layout.

## Out of Scope

- interactive filtering, zooming, or sample selection;
- external plotting libraries or network-loaded assets;
- changing peak calling, FRiP, TSS, or motif statistical methods;
- replacing the existing raw QC files;
- pass/fail classification or universal biological thresholds.
