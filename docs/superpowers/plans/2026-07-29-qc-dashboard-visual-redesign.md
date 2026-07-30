# QC Dashboard Visual Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the consolidated offline QC report into a publication-style dashboard with explicit units, value labels, target-consistent colors, normalized 250-bp distributions, peak-coverage ECDFs, direct sample labels, and a cognate-aware motif heatmap.

**Architecture:** Keep input parsing, validation, joining, and output contracts in `bin/qc_dashboard.py`. Add `bin/qc_dashboard_visuals.py` for deterministic data-to-geometry transforms and dependency-free inline SVG rendering, then pass the existing per-peak fragment-count output through the QC subworkflow into the dashboard. Preserve full-resolution JSON/TSV data while rendering compact, collapsed HTML tables.

**Tech Stack:** Python 3.12 standard library, Nextflow DSL2, inline HTML/CSS/SVG, `unittest`, shell/Python integration tests.

## Global Constraints

- Keep the dashboard a deterministic, self-contained HTML file with no external scripts, styles, fonts, data, JavaScript dependency, or plotting library.
- Keep the report descriptive; do not add biological pass/fail classifications or universal QC thresholds.
- Use shared, left-closed 250-bp bins for both insert size and peak width.
- Normalize insert-size bins to percent of observed read pairs and peak-width bins to percent of called peaks within each sample.
- Render fragments per peak as an ECDF with a labeled zero origin and a logarithmic positive domain.
- Use IgG grey, CTCF blue, GATA1 green, and RUNX1 orange consistently; assign future targets deterministic accessible colors.
- Supplement color with markers, dashes, hatches, visible labels, and accessible SVG titles.
- Format visible HTML table metrics to three significant digits while preserving full precision in TSV and JSON.
- Retain full-resolution insert-size and peak-width arrays in `qc_summary.json`.
- Keep `top_motifs.tsv` limited to ten motifs per target while using complete AME tables to construct the HTML heatmap.
- Increment the dashboard generator version from `1.0.0` to `1.1.0`; keep schema version `1`.
- Missing optional inputs render an unavailable state; present malformed inputs fail with a filename and contract reason.
- Publish the same five dashboard files atomically under `<outdir>/reports/qc_dashboard/`.

---

## File Structure

- Create `bin/qc_dashboard_visuals.py`: pure formatting, binning, ECDF, color/style, label-packing, and SVG-rendering helpers.
- Create `tests/unit/test_qc_dashboard_visuals.py`: direct tests for all visual transforms and SVG contracts.
- Modify `bin/qc_dashboard.py`: parse fragments-per-peak data and complete AME data, extend the internal/public model, assemble redesigned sections, collapse tables, and update generator metadata.
- Modify `tests/unit/test_qc_dashboard.py`: test new parsers, joins, public JSON, heatmap selection, report sections, compact tables, and CLI behavior.
- Modify `subworkflows/local/qc.nf`: collect `*.fragments_per_peak.tsv` from `PEAK_QC` and pass it to `QC_DASHBOARD`.
- Modify `modules/local/qc_dashboard.nf`: stage fragments-per-peak inputs and update the version record.
- Modify `tests/integration/test_qc.sh`: verify Nextflow wiring and directly exercise the staged fragments-per-peak contract.
- Modify `tests/integration/test_e2e.sh`: require redesigned report semantics and compact outputs in the offline workflow.
- Modify `tests/data/e2e/fakebin/fake_bio_tool.py`: emit complete multi-motif AME fixture rows needed by the heatmap.
- Create `tests/render_qc_dashboard_fixture.py`: generate a representative 12-sample HTML report for repeatable visual QA.
- Modify `README.md`: document the redesigned sections, 250-bp transformations, fragments-per-peak metric, motif heatmap, table precision, and generator version.
- Modify `CHANGELOG.md`: record the dashboard visual redesign.
- Modify `tests/integration/test_docs.sh`: assert the new documented generator version and report semantics.

### Task 1: Pure visual transforms and formatting

**Files:**
- Create: `bin/qc_dashboard_visuals.py`
- Create: `tests/unit/test_qc_dashboard_visuals.py`

**Interfaces:**
- Produces:
  - `TARGET_COLORS: Mapping[str, str]`
  - `target_color(assay_target: object, *, is_control: bool = False) -> str`
  - `SeriesStyle` dataclass with `color`, `dash`, `marker`, and `is_control`
  - `series_style(sample_id: str, assay_target: object, *, is_control: bool) -> SeriesStyle`
  - `format_significant(value: object, *, exact: bool = False, compact: bool = True) -> str`
  - `bin_weighted_series(series: Mapping[str, Sequence[tuple[int, int]]], *, bin_size: int = 250) -> dict[str, list[dict[str, float | int]]]`
  - `histogram_ecdf(histogram: Sequence[tuple[int, int]]) -> list[dict[str, float | int]]`
  - `pack_endpoint_labels(endpoints: Sequence[tuple[str, float]], *, lower: float, upper: float, minimum_gap: float) -> dict[str, float]`

- [ ] **Step 1: Write failing color and three-significant-digit tests**

Add `unittest.TestCase` tests:

```python
def test_target_colors_are_stable_and_controls_are_grey(self):
    self.assertEqual(visuals.target_color("IgG", is_control=True), "#8F8D87")
    self.assertEqual(visuals.target_color("CTCF"), "#2F78D1")
    self.assertEqual(visuals.target_color("GATA1"), "#1FAE7A")
    self.assertEqual(visuals.target_color("RUNX1"), "#F06432")
    self.assertEqual(visuals.target_color("NEW_TF"), visuals.target_color("NEW_TF"))

def test_html_numbers_use_three_significant_digits_without_rounding_exact_fields(self):
    self.assertEqual(visuals.format_significant(17_432_100), "17.4M")
    self.assertEqual(visuals.format_significant(0.012345), "0.0123")
    self.assertEqual(visuals.format_significant(0.000012345), "1.23e-05")
    self.assertEqual(visuals.format_significant(1250, exact=True), "1250")
```

- [ ] **Step 2: Run the tests and verify RED**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard_visuals.VisualFormattingTests -v
```

Expected: FAIL because `qc_dashboard_visuals` does not exist.

- [ ] **Step 3: Implement deterministic colors, styles, and formatting**

Create the module with:

```python
TARGET_COLORS = {
    "IGG": "#8F8D87",
    "CTCF": "#2F78D1",
    "GATA1": "#1FAE7A",
    "RUNX1": "#F06432",
}
FALLBACK_COLORS = (
    "#7C3AED", "#0F766E", "#BE123C", "#A16207",
    "#0369A1", "#6D28D9", "#047857", "#B45309",
)

@dataclass(frozen=True)
class SeriesStyle:
    color: str
    dash: str
    marker: str
    is_control: bool
```

Use `hashlib.sha256(normalized_target.encode()).digest()[0]` to select a
fallback color; never use Python's randomized `hash()`. Derive dash/marker
from a SHA-256 digest of `sample_id`, using a square plus dashed stroke for
controls and a circle plus deterministic dash for targets.

`format_significant` must return `NA` for `None`, lowercase booleans, exact
strings/integers when requested, `K/M/B` compact values at powers of 1000,
three significant digits for finite metrics, and scientific notation below
`1e-4`.

- [ ] **Step 4: Write failing 250-bp binning and ECDF tests**

```python
def test_shared_250bp_bins_include_zero_bins_and_normalize_each_sample(self):
    result = visuals.bin_weighted_series({
        "A": [(0, 1), (249, 1), (250, 2), (501, 1)],
        "B": [(260, 4)],
    })
    self.assertEqual(
        [(row["bin_start"], row["bin_end"]) for row in result["A"]],
        [(0, 250), (250, 500), (500, 750)],
    )
    self.assertEqual(
        [row["count"] for row in result["B"]],
        [0, 4, 0],
    )
    self.assertAlmostEqual(sum(row["percent"] for row in result["A"]), 100.0)

def test_fragment_histogram_ecdf_retains_zero_and_reaches_one_hundred_percent(self):
    rows = visuals.histogram_ecdf([(0, 2), (1, 1), (10, 1)])
    self.assertEqual(rows[0], {"value": 0, "count": 2, "cumulative_percent": 50.0})
    self.assertEqual(rows[-1]["cumulative_percent"], 100.0)
```

- [ ] **Step 5: Run the transform tests and verify RED**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard_visuals.DistributionTransformTests -v
```

Expected: FAIL because the transform functions are missing.

- [ ] **Step 6: Implement shared bins, ECDFs, and label packing**

`bin_weighted_series` must validate positive `bin_size`, non-negative integer
positions/counts, create bins from zero through the largest observed value,
fill absent sample bins with zero, and calculate `percent = 100 * count /
sample_total` or `0.0` for an empty sample.

`histogram_ecdf` must combine duplicate values, reject negative values/counts,
omit zero-count rows, and return cumulative percentages in numeric order.

`pack_endpoint_labels` must sort by requested y coordinate, perform a forward
minimum-gap pass, shift overflow upward, perform a reverse pass, and constrain
all labels to `[lower, upper]`. Reject an impossible layout when
`minimum_gap * (n - 1) > upper - lower`.

- [ ] **Step 7: Run Task 1 tests and commit**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard_visuals -v
git diff --check
```

Expected: all Task 1 tests pass.

Commit:

```bash
git add bin/qc_dashboard_visuals.py tests/unit/test_qc_dashboard_visuals.py
git commit -m "feat: add dashboard visual transforms"
```

### Task 2: Fragments-per-peak parser, model, and Nextflow transport

**Files:**
- Modify: `bin/qc_dashboard.py`
- Modify: `tests/unit/test_qc_dashboard.py`
- Modify: `subworkflows/local/qc.nf`
- Modify: `modules/local/qc_dashboard.nf`
- Modify: `tests/integration/test_qc.sh`

**Interfaces:**
- Produces:
  - `read_fragments_per_peak_distributions(paths: Sequence[Path]) -> dict[str, list[dict[str, int]]]`
  - internal `sample["peak"]["fragments_per_peak_distribution"]`
  - public JSON `sample["peak"]["fragments_per_peak_distribution"]`
- Consumes the exact producer header:
  - `chrom start end peak_name width score signal_value fragment_count`

- [ ] **Step 1: Write failing parser and public-contract tests**

Add tests that create:

```text
chrom	start	end	peak_name	width	score	signal_value	fragment_count
chr1	0	250	p1	250	10	3.5	0
chr1	500	1000	p2	500	12	5.0	4
chr2	0	750	p3	750	8	2.0	4
```

Require:

```python
self.assertEqual(
    qc.read_fragments_per_peak_distributions([path]),
    {"S1": [
        {"fragment_count": 0, "peak_count": 1},
        {"fragment_count": 4, "peak_count": 2},
    ]},
)
```

Also require rejection of a wrong header, negative/fractional counts, duplicate
sample files, and a filename that does not end in
`.peak_qc.fragments_per_peak.tsv`. Extend the report-data test so the public
JSON preserves the compact histogram and does not contain chromosome or peak
name rows.

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard.MetadataAndTableParserTests \
  tests.unit.test_qc_dashboard.DashboardOutputTests -v
```

Expected: FAIL because the parser/model field does not exist.

- [ ] **Step 3: Implement strict parsing and model joins**

Derive sample ID with the suffix
`.peak_qc.fragments_per_peak.tsv`. Require the exact eight-column header and
parse only `fragment_count`; aggregate with `collections.Counter`. A
header-only file is a valid empty target.

Add optional keyword input to `build_report_data`:

```python
fragments_per_peak: Mapping[str, Sequence[Mapping[str, object]]] | None = None
```

Reject unknown sample IDs. For each target with peak metrics, set the histogram
and distinguish computed versus empty. Controls retain `None` and
`not_applicable`.

Add to `_public_sample`:

```python
public_peak["fragments_per_peak_distribution"] = [
    {
        "fragment_count": row.get("fragment_count"),
        "peak_count": row.get("peak_count"),
    }
    for row in distribution
]
```

- [ ] **Step 4: Write failing Nextflow wiring checks**

Extend `tests/integration/test_qc.sh` to require:

```python
"QC collects fragments-per-peak outputs":
    ".map { meta, json, tsv, histogram, perPeak -> perPeak }" in qc_subworkflow,
"dashboard stages fragments-per-peak files":
    "path peak_fragment_counts" in dashboard
    and "dashboard_inputs/peak/fragments??/*" in dashboard,
```

Extend the extracted dashboard preprocess fixture with two staged
`*.fragments_per_peak.tsv` files and require successful identity validation.

- [ ] **Step 5: Run integration checks and verify RED**

Run:

```bash
bash tests/integration/test_qc.sh
```

Expected: FAIL on the new transport assertions.

- [ ] **Step 6: Implement Nextflow collection and staging**

In `subworkflows/local/qc.nf`, add:

```groovy
peak_fragment_count_files = PEAK_QC.out.qc
    .map { meta, json, tsv, histogram, perPeak -> perPeak }
    .reduce([files: []]) { holder, path ->
        [files: holder.files + [path]]
    }
    .map { holder -> holder.files }
```

Pass it immediately after `peak_width_histogram_files` to `QC_DASHBOARD`.

In `modules/local/qc_dashboard.nf`, add:

```groovy
path peak_fragment_counts,
    stageAs: 'dashboard_inputs/peak/fragments??/*'
```

The CLI continues to use `--peak-dir`; `main()` finds
`*.peak_qc.fragments_per_peak.tsv` recursively and passes the parsed mapping
to `build_report_data`.

- [ ] **Step 7: Run Task 2 checks and commit**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard -q
bash tests/integration/test_qc.sh
git diff --check
```

Commit:

```bash
git add bin/qc_dashboard.py tests/unit/test_qc_dashboard.py \
  subworkflows/local/qc.nf modules/local/qc_dashboard.nf \
  tests/integration/test_qc.sh
git commit -m "feat: carry peak coverage into dashboard"
```

### Task 3: Complete AME data and cognate heatmap matrix

**Files:**
- Modify: `bin/qc_dashboard.py`
- Modify: `tests/unit/test_qc_dashboard.py`

**Interfaces:**
- Produces:
  - `read_all_ame(path: Path) -> list[dict[str, object]]`
  - `motif_tokens(*values: object) -> frozenset[str]`
  - `is_cognate_motif(expected_motif: object, motif_id: object, motif_alt_id: object) -> bool`
  - `build_motif_heatmap(samples: Sequence[Mapping[str, object]], *, noncognate_limit: int = 15, significance_cap: float = 60.0) -> dict[str, object]`
- Preserves `read_top_ame(..., limit=10)` and `top_motifs.tsv`.

- [ ] **Step 1: Write failing full-AME and token-aware matching tests**

Require `read_all_ame` to retain more than ten rows while `read_top_ame`
remains limited. Add:

```python
self.assertTrue(qc.is_cognate_motif("GATA1", "MA0140.2", "GATA1::TAL1"))
self.assertTrue(qc.is_cognate_motif("GATA1", "TAL1::GATA1", "complex"))
self.assertFalse(qc.is_cognate_motif("GATA1", "MA9999", "GATA10"))
```

Use punctuation-delimited, case-insensitive tokens from both ID and alternate
name.

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard.TssAndAmeTests -v
```

Expected: FAIL because full AME and cognate helpers are absent.

- [ ] **Step 3: Implement full AME parsing without changing the top-ten export**

Refactor the existing serializer:

```python
def _ame_record_dict(record):
    return {
        "motif_id": record.motif_id,
        "motif_alt_id": record.motif_alt_id,
        "adjusted_p_value": record.adjusted_p_value,
        "p_value": record.p_value,
        "effect": record.effect,
        "positive_sequences": record.positive_sequences,
        "rank": record.rank,
    }
```

`read_all_ame` returns every non-`__NO_PEAKS__` record in deterministic
significance/rank/identity order. `read_top_ame` slices that result. Change the
directory loader to return both complete and top-ten mappings without reading
each AME file twice. Store complete records internally as
`sample["ame_motifs"]`; keep public `top_motifs` bounded.

- [ ] **Step 4: Write failing cohort-selection and heatmap-cell tests**

Create samples with:

- 20 noncognate motifs;
- `GATA1`, `GATA1::TAL1`, and `GATA10`;
- adjusted values `0`, `0.001`, `0.2`, and missing;
- two target samples and one IgG control.

Require:

```python
matrix = qc.build_motif_heatmap(samples)
self.assertNotIn("IgG_sample", matrix["sample_ids"])
self.assertIn(("MA_GATA1", "GATA1::TAL1"), matrix["motif_keys"])
self.assertNotIn(("MA_GATA10", "GATA10"), matrix["forced_cognate_keys"])
self.assertEqual(len(matrix["noncognate_keys"]), 15)
self.assertEqual(matrix["cells"][("GATA_sample", ("MA_ZERO", "ZERO"))]["score"], 60.0)
self.assertEqual(matrix["cells"][("GATA_sample", ("MA_NS", "NS"))]["label"], "ns")
```

- [ ] **Step 5: Implement deterministic motif selection and cells**

Tokenize with:

```python
re.findall(r"[A-Za-z0-9]+", text.casefold())
```

An expected TF is cognate when its complete normalized token occurs in either
motif field. Select all cognate motif keys first, grouped by expected TF, then
the 15 remaining keys with the smallest adjusted value observed in any sample.
Tie-break by alternate name then motif ID.

For each sample/key cell:

- `adjusted == 0`: score `60.0`, label `>60`;
- `0 < adjusted <= 0.1`: score `min(-log10(adjusted), 60)`;
- `adjusted > 0.1` or absent: label `ns`;
- `cognate`: `outlined=True`.

- [ ] **Step 6: Run Task 3 tests and commit**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard.TssAndAmeTests -v
git diff --check
```

Commit:

```bash
git add bin/qc_dashboard.py tests/unit/test_qc_dashboard.py
git commit -m "feat: build cognate motif heatmap data"
```

### Task 4: SVG primitives and sequencing/peak small multiples

**Files:**
- Modify: `bin/qc_dashboard_visuals.py`
- Modify: `tests/unit/test_qc_dashboard_visuals.py`
- Modify: `bin/qc_dashboard.py`
- Modify: `tests/unit/test_qc_dashboard.py`

**Interfaces:**
- Produces:
  - `BarMetric` dataclass describing key, title, axis label, scale, and formatter
  - `render_bar_panel(title: str, rows: Sequence[Mapping[str, object]], *, value_key: str, axis_label: str, value_multiplier: float = 1.0, log10_axis: bool = False) -> str`
  - `render_range_panel(...) -> str`
  - `render_scatter_panel(...) -> str`
  - `render_panel_grid(panels: Sequence[str], *, aria_label: str) -> str`
- Consumes rows with `sample_id`, `assay_target`, `is_control`, and numeric values.

- [ ] **Step 1: Write failing SVG accessibility, axis, and value-label tests**

Require one bar panel to contain:

```python
self.assertIn('aria-label="Mapped reads (%)"', panel)
self.assertIn('class="axis-title axis-title-y"', panel)
self.assertIn('class="bar-value"', panel)
self.assertIn(">82.9%<", panel)
self.assertIn('data-assay-target="CTCF"', panel)
self.assertIn('fill="#2F78D1"', panel)
```

Add tests for a log10 usable-fragment axis with labeled powers, an IQR/range
panel, a labeled scatterplot, empty panels, long sample labels, and 100-sample
horizontal scrolling.

- [ ] **Step 2: Run focused visual tests and verify RED**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard_visuals.BarAndScatterRenderingTests -v
```

Expected: FAIL because the renderers are absent.

- [ ] **Step 3: Implement reusable axes and bar geometry**

Render at least four y ticks plus baseline, all with visible tick labels. Bar
height uses either a linear domain beginning at zero or explicit log10 values;
the axis title includes `(log10)` for transformed metrics. Use
`target_color/series_style`, visible bar-value text, exact values in `<title>`,
and rotated sample labels only when required by measured text length.

`render_range_panel` draws minimum–maximum lines, Q25–Q75 thick segments, and
median circles. `render_scatter_panel` draws one point per sample plus a direct
label and leader line; both axes receive titles and ticks.

- [ ] **Step 4: Write failing report-section tests for all approved panels**

In `tests/unit/test_qc_dashboard.py`, require these titles:

```text
Assigned read pairs
Barcode balance within library
Mapped reads
Usable fragments after filtering
PCR duplication
End-to-end usable yield
Peak count
Fraction of reads in peaks
Total bases covered by peaks
Peak width median and range
Peak count vs usable fragments
TSS enrichment score
```

Require the demultiplex denominator to be assigned reads and end-to-end yield
to equal `100 * mapq_filtered_fragments / sample_assigned_reads`, with `NA`
when the denominator is zero or missing.

- [ ] **Step 5: Run report tests and verify RED**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard.DashboardOutputTests -v
```

Expected: FAIL because `render_dashboard` still renders three legacy bars.

- [ ] **Step 6: Assemble sequencing and peak panel grids**

Add a pure derived-row helper in `qc_dashboard.py` that computes:

```python
assigned_read_pairs_millions = sample_assigned_reads / 1_000_000
barcode_balance_percent = 100 * sample_assignment_fraction
usable_fragments_millions = mapq_filtered_fragments / 1_000_000
end_to_end_yield_percent = 100 * mapq_filtered_fragments / sample_assigned_reads
peak_count_thousands = peak_count / 1_000
frip_percent = 100 * frip
covered_megabases = total_covered_bases / 1_000_000
```

Guard every division and retain `None` for unavailable values. Assemble
responsive `.panel-grid` sections using the SVG helpers. Keep controls in
technical panels and exclude controls from biological peak panels.

- [ ] **Step 7: Run Task 4 tests and commit**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard_visuals \
  tests.unit.test_qc_dashboard -q
git diff --check
```

Commit:

```bash
git add bin/qc_dashboard_visuals.py tests/unit/test_qc_dashboard_visuals.py \
  bin/qc_dashboard.py tests/unit/test_qc_dashboard.py
git commit -m "feat: render dashboard QC panel grids"
```

### Task 5: Target-colored distributions, ECDFs, and direct labels

**Files:**
- Modify: `bin/qc_dashboard_visuals.py`
- Modify: `tests/unit/test_qc_dashboard_visuals.py`
- Modify: `bin/qc_dashboard.py`
- Modify: `tests/unit/test_qc_dashboard.py`

**Interfaces:**
- Produces:
  - `render_binned_distribution(title: str, series: Mapping[str, Sequence[Mapping[str, object]]], metadata: Mapping[str, Mapping[str, object]], *, x_axis_label: str) -> str`
  - `render_ecdf(title: str, series: Mapping[str, Sequence[Mapping[str, object]]], metadata: Mapping[str, Mapping[str, object]], *, x_axis_label: str, zero_origin: bool) -> str`
  - `render_profile_chart(...) -> str`
- Uses `pack_endpoint_labels` and `series_style` from Task 1.

- [ ] **Step 1: Write failing distribution geometry and label tests**

Require binned series to:

- use `percent` on the y-axis;
- place x at bin midpoint;
- use one target color for all samples of the same TF;
- retain distinct dashes/markers;
- emit one `class="endpoint-label"` and optional leader per sample;
- keep packed label y coordinates inside the plot with the minimum gap.

Require the fragments-per-peak ECDF to contain a separate `0` x tick and
positive log ticks `1`, `10`, `100`, and so on.

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard_visuals.DistributionRenderingTests -v
```

Expected: FAIL because the new renderers are absent.

- [ ] **Step 3: Implement normalized distribution and ECDF renderers**

`render_binned_distribution` uses the shared bins returned by Task 1 and
renders percent traces on a linear x-axis. `render_ecdf` maps zero to a
separate origin slot and maps positive counts with `log10(value)`. Both
renderers call `pack_endpoint_labels`, draw leader lines when displaced, add
exact-value SVG titles, and retain a compact legend.

If endpoint packing is impossible at the fixed plot height, increase the SVG
height deterministically instead of dropping labels.

- [ ] **Step 4: Write failing dashboard distribution and TSS tests**

Require:

- insert and peak-width plots/table rows use `[0,250)`, `[250,500)`, ...;
- all samples share bins through the cohort maximum;
- binned percentages sum to 100 per nonempty sample;
- fragments-per-peak uses an ECDF and reports zero-peak percentage;
- TSS profiles use target colors, a visible x=0 reference, RPKM y-axis, and
  direct sample labels.

- [ ] **Step 5: Run report tests and verify RED**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard.DashboardOutputTests -v
```

Expected: FAIL because the report still renders raw-count distributions and
legacy TSS styles.

- [ ] **Step 6: Replace legacy distribution/TSS rendering**

Build weighted series:

```python
insert_series[sample_id] = [
    (row["insert_size"], row["pair_count"]) for row in distribution
]
peak_width_series[sample_id] = [
    (row["width"], row["peak_count"]) for row in distribution
]
```

Pass them to `bin_weighted_series(bin_size=250)`. Convert the compact
fragments-per-peak histograms with `histogram_ecdf`. Add `assay_target` to TSS
profile metadata and render profiles through `render_profile_chart`.

- [ ] **Step 7: Run Task 5 tests and commit**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard_visuals \
  tests.unit.test_qc_dashboard -q
git diff --check
```

Commit:

```bash
git add bin/qc_dashboard_visuals.py tests/unit/test_qc_dashboard_visuals.py \
  bin/qc_dashboard.py tests/unit/test_qc_dashboard.py
git commit -m "feat: render normalized dashboard distributions"
```

### Task 6: Cognate-aware motif heatmap rendering

**Files:**
- Modify: `bin/qc_dashboard_visuals.py`
- Modify: `tests/unit/test_qc_dashboard_visuals.py`
- Modify: `bin/qc_dashboard.py`
- Modify: `tests/unit/test_qc_dashboard.py`

**Interfaces:**
- Produces:
  - `render_motif_heatmap(matrix: Mapping[str, object]) -> str`
- Consumes the exact `build_motif_heatmap` result from Task 3.

- [ ] **Step 1: Write failing heatmap rendering tests**

Require:

```python
self.assertIn('aria-label="Motif enrichment heatmap"', heatmap)
self.assertIn('class="heatmap-cell cognate"', heatmap)
self.assertIn(">60<", heatmap)
self.assertIn(">ns<", heatmap)
self.assertIn("−log10 adjusted p-value (capped at 60)", heatmap)
self.assertIn('class="heatmap-scroll"', heatmap)
```

Also require target-colored column labels, exact adjusted values in cell
titles, unique motif/sample accessible names, and a minimum cell width.

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard_visuals.MotifHeatmapRenderingTests -v
```

Expected: FAIL because the renderer is absent.

- [ ] **Step 3: Implement the static heatmap**

Use a white-to-blue sequential scale computed from score/60. Draw a dark
outline for cognate cells. Render row labels, rotated target-colored sample
labels, centered cell text, accessible `<title>` elements, and an inline
0–60 legend. Set SVG width from `row_label_width + cell_width * sample_count`
and place it inside `.heatmap-scroll`.

- [ ] **Step 4: Write failing dashboard heatmap/table tests**

Require the motif section to contain the heatmap and a compact table with only
displayed matrix cells. Confirm that an eleventh noncognate top motif can enter
the heatmap while `top_motifs.tsv` remains limited to ten. Confirm controls do
not become columns.

- [ ] **Step 5: Assemble the motif section**

Call `build_motif_heatmap(_samples(data))`, render it, and construct compact
rows with:

```text
sample_id, motif_id, motif_alt_id, rank, adjusted_p_value,
transformed_significance, cognate
```

Do not place `ame_motifs` into public JSON. Continue to serialize only
`top_motifs`.

- [ ] **Step 6: Run Task 6 tests and commit**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard_visuals \
  tests.unit.test_qc_dashboard -q
git diff --check
```

Commit:

```bash
git add bin/qc_dashboard_visuals.py tests/unit/test_qc_dashboard_visuals.py \
  bin/qc_dashboard.py tests/unit/test_qc_dashboard.py
git commit -m "feat: render cognate motif heatmap"
```

### Task 7: Collapsed compact tables, versioning, and responsive layout

**Files:**
- Modify: `bin/qc_dashboard.py`
- Modify: `tests/unit/test_qc_dashboard.py`
- Modify: `modules/local/qc_dashboard.nf`

**Interfaces:**
- Produces:
  - `render_table(..., exact_keys: Collection[str] = ()) -> str`
  - `render_data_details(title: str, content: str, *, open_by_default: bool = False) -> str`
- Uses `format_significant` from Task 1.

- [ ] **Step 1: Write failing precision and collapsed-table tests**

Require:

```python
self.assertIn(">17.4M<", dashboard)
self.assertIn(">0.0123<", dashboard)
self.assertIn("<details", dashboard)
self.assertIn("<summary>View binned insert-size data</summary>", dashboard)
self.assertNotIn("<details open", dashboard)
```

Assert sample IDs, ranks, and `250` bin boundaries remain exact. Serialize the
same data and assert JSON/TSV retain `17432100` and `0.012345` without HTML
rounding.

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tests.unit.test_qc_dashboard.DashboardOutputTests -v
```

Expected: FAIL because current tables use six significant digits and remain
expanded.

- [ ] **Step 3: Implement presentation-only table formatting**

`render_table` calls:

```python
format_significant(row.get(key), exact=key in exact_keys)
```

Pass exact keys for sample/library IDs, motif IDs, categories, ranks, and bin
boundaries. Wrap detailed metric, distribution, TSS, and motif tables in
`render_data_details`; leave run overview, input availability, and warnings
visible.

- [ ] **Step 4: Write failing version/layout tests**

Require:

```python
self.assertEqual(qc.GENERATOR_VERSION, "1.1.0")
self.assertIn("qc_dashboard.py: 1.1.0", dashboard_module)
self.assertIn(".panel-grid{", html)
self.assertIn("@media (max-width:", html)
self.assertIn("@media print{", html)
```

Keep `SCHEMA_VERSION == 1` and the five output filenames unchanged.

- [ ] **Step 5: Update version metadata and responsive CSS**

Set `GENERATOR_VERSION = "1.1.0"` and the Nextflow YAML record to `1.1.0`.
Add responsive two/three-column panel grids, narrow-width single-column rules,
print rules that avoid section splitting, chart/heatmap scroll containers, and
visible focus states for scrollable/table detail regions.

- [ ] **Step 6: Run Task 7 tests and commit**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s tests/unit -p 'test_*.py' -q
bash -n tests/integration/test_qc.sh
git diff --check
```

Commit:

```bash
git add bin/qc_dashboard.py tests/unit/test_qc_dashboard.py \
  modules/local/qc_dashboard.nf
git commit -m "feat: compact dashboard data tables"
```

### Task 8: Offline end-to-end coverage and documentation

**Files:**
- Modify: `tests/integration/test_e2e.sh`
- Modify: `tests/data/e2e/fakebin/fake_bio_tool.py`
- Modify: `tests/integration/test_docs.sh`
- Modify: `README.md`
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes the five unchanged dashboard output paths.
- Documents generator `1.1.0`, 250-bp bins, coverage ECDF, target palette,
  motif heatmap, and three-significant-digit HTML presentation.

- [ ] **Step 1: Write failing end-to-end assertions and richer fake AME output**

Update fake AME output to include more than 15 motifs, including:

```text
MA0139.1	CTCF
MA1929.1	CTCF::ZNF143
MA0140.2	GATA1::TAL1
MA1356.1	TAL1::GATA1
MA0002.2	RUNX1
MA9999.1	GATA10
```

Mark the fixture as complete database output. Extend E2E assertions for:

- all small-multiple titles;
- `250 bp bins`;
- `Fragments per peak`;
- `Motif enrichment heatmap`;
- `class="endpoint-label"`;
- collapsed `<details>`;
- no external `http://` or `https://`;
- unchanged raw insert/width JSON values;
- compact fragments-per-peak JSON histogram.

- [ ] **Step 2: Run E2E checks and verify RED**

Run:

```bash
bash tests/integration/test_e2e.sh
```

Expected: FAIL on the new fixture/report assertions.

- [ ] **Step 3: Complete E2E fixtures and assertions**

Update fake output and expected report checks without weakening existing
resume, output-path, MultiQC, or atomic-publication assertions. Keep fixture
sizes small enough for local execution.

- [ ] **Step 4: Update README, changelog, and documentation checks**

Document:

- the six sequencing panels and five peak panels;
- normalized shared 250-bp insert/width bins;
- fragments-per-peak source and zero handling;
- target colors and direct labels;
- complete-AME top-15-plus-cognate heatmap selection;
- three-significant-digit HTML display versus full-precision TSV/JSON;
- generator version `1.1.0`, schema version `1`.

Update `tests/integration/test_docs.sh` to require those exact contract phrases.

- [ ] **Step 5: Run Task 8 checks and commit**

Run:

```bash
bash tests/integration/test_docs.sh
bash tests/integration/test_qc.sh
bash tests/integration/test_e2e.sh
git diff --check
```

Commit:

```bash
git add tests/integration/test_e2e.sh \
  tests/data/e2e/fakebin/fake_bio_tool.py \
  tests/integration/test_docs.sh README.md CHANGELOG.md
git commit -m "docs: describe redesigned QC dashboard"
```

### Task 9: Representative fixture, visual QA, and final verification

**Files:**
- Create: `tests/render_qc_dashboard_fixture.py`
- Test: all dashboard and workflow tests

**Interfaces:**
- Produces a local QA artifact at a caller-supplied output path.
- Uses public `build_report_data`/`render_dashboard` interfaces; does not add a
  production CLI mode.

- [ ] **Step 1: Create a deterministic 12-sample rendering fixture**

Build three IgG controls, three CTCF targets, three GATA1 targets, and three
RUNX1 targets with:

- distinct read depth, mapping, duplicate, FRiP, peak count, width, and TSS
  values;
- insert/width values crossing several 250-bp boundaries;
- zero and positive fragments-per-peak counts;
- overlapping line endpoints that exercise label packing;
- significant, nonsignificant, zero-p-value, exact, and complex cognate motifs.

Expose:

```bash
python3 tests/render_qc_dashboard_fixture.py /tmp/qc-dashboard-visual.html
```

- [ ] **Step 2: Add a smoke test for the fixture**

Add a test in `tests/unit/test_qc_dashboard.py` that invokes the fixture builder
and requires a nonempty self-contained HTML document with 12 samples and all
approved panels.

- [ ] **Step 3: Run the complete automated suite**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s tests/unit -p 'test_*.py' -q
bash tests/integration/test_manifest.sh
bash tests/integration/test_demultiplex.sh
bash tests/integration/test_alignment.sh
bash tests/integration/test_peaks.sh
bash tests/integration/test_qc.sh
bash tests/integration/test_motifs.sh
bash tests/integration/test_docs.sh
bash tests/integration/test_e2e.sh
git diff --check
```

Expected: all unit and integration checks pass. If Nextflow is unavailable,
integration scripts may report their existing explicit DSL2 skip while all
static/direct checks must pass.

- [ ] **Step 4: Render and inspect desktop, narrow, and print layouts**

Generate:

```bash
python3 tests/render_qc_dashboard_fixture.py /tmp/qc-dashboard-visual.html
```

Open the local file in a browser and inspect at:

- 1440 × 1000 desktop;
- 390 × 844 narrow viewport;
- print preview.

Verify:

- every bar has a readable value and y-axis title;
- target colors remain consistent across every section;
- direct labels do not overlap and leader lines point to the correct series;
- 250-bp binned lines and compact tables agree;
- fragments-per-peak shows zero and positive log ticks;
- the motif heatmap scrolls without expanding page width;
- details are collapsed and keyboard-focusable;
- no chart, table, or section is clipped in print preview.

Record screenshots under `/tmp` only; do not commit generated images.

- [ ] **Step 5: Review the complete diff against the specification**

Run:

```bash
git diff --stat origin/agent/bulk-nanocut-tag-pipeline...HEAD
git diff --check origin/agent/bulk-nanocut-tag-pipeline...HEAD
git status -sb
```

Confirm every requirement in
`docs/superpowers/specs/2026-07-29-qc-dashboard-visual-redesign-design.md`
maps to an implemented test and no unrelated files are changed.

- [ ] **Step 6: Commit the visual fixture**

```bash
git add tests/render_qc_dashboard_fixture.py tests/unit/test_qc_dashboard.py
git commit -m "test: add dashboard visual QA fixture"
```

- [ ] **Step 7: Push and update draft PR #11**

```bash
git push origin codex/qc-dashboard-visual-redesign
gh pr view 11 --json state,isDraft,url,headRefOid
```

Keep PR #11 in draft until the full automated suite and visual QA pass.
