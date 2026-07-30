# Consolidated QC Dashboard Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a self-contained run-level QC dashboard and stable TSV/JSON summaries covering demultiplexing, alignment, peaks, FRiP, TSS enrichment, and motif enrichment while preserving MultiQC.

**Architecture:** A repository-owned Python program will parse the existing QC artifacts, join them to manifest-derived sample metadata, calculate a documented TSS score, and render deterministic TSV, JSON, and self-contained HTML outputs. A new Nextflow module will collect those artifacts after the existing QC processes, publish the dashboard, and expose it through the top-level workflow and completion summary.

**Tech Stack:** Python 3.12 standard library, Nextflow DSL2, deepTools 3.5.5, inline HTML/CSS/SVG, shell/Python integration tests.

## Global Constraints

- Preserve `<outdir>/reports/multiqc/multiqc_report.html` and `<outdir>/reports/summary/combined_target_qc.tsv`.
- Publish the new outputs under `<outdir>/reports/qc_dashboard/`.
- The HTML report must be self-contained and must not load external scripts, styles, fonts, or data.
- The first release is descriptive: do not assign biological pass/fail labels or fixed QC thresholds.
- Include all libraries in technical QC and visually distinguish IgG controls from target libraries.
- Exclude IgG controls from expected-motif interpretation while retaining their technical QC.
- Missing optional metrics become blank TSV values and JSON `null`, with a visible warning; they do not become zero.
- Duplicate sample identifiers, contradictory metadata, and incompatible input schemas remain fatal.
- Calculate TSS enrichment as the center-bin signal divided by the mean signal in the terminal 100 bp at both flanks.
- Report the expected motif and up to ten top-ranked AME motifs for every target sample.
- Use only the existing `envs/python.yml` runtime for the dashboard; add no plotting or web dependencies.

---

## File Structure

- Create `bin/qc_dashboard.py`: parse QC artifacts, join sample records, calculate TSS enrichment, write deterministic TSV/JSON files, and render inline SVG/HTML.
- Create `tests/unit/test_qc_dashboard.py`: direct tests of parsing, validation, calculations, serialization, and self-contained rendering.
- Create `modules/local/qc_dashboard.nf`: stage collected QC artifacts, invoke the Python generator, publish five dashboard files, and emit a version record.
- Modify `modules/local/tss_enrichment.nf`: ask `plotProfile` for its average-profile data table and emit it with the existing TSS outputs.
- Modify `subworkflows/local/motifs.nf`: preserve the existing AME result/status outputs for downstream dashboard consumption.
- Modify `subworkflows/local/qc.nf`: create validated sample metadata JSON, collect detailed QC inputs, invoke `QC_DASHBOARD`, and expose its outputs.
- Modify `main.nf`: pass AME inputs to QC, expose the dashboard outputs, and add the dashboard to `run_summary.txt`.
- Modify `tests/integration/test_qc.sh`: assert the expanded QC contracts and directly exercise dashboard generation with missing optional data.
- Modify `tests/integration/test_e2e.sh`: assert the five published dashboard files, their content, and resume behavior.
- Create `tests/data/e2e/tss.bed`: deterministic BED6 TSS fixture.
- Create `tests/data/e2e/fakebin/computeMatrix` and `tests/data/e2e/fakebin/plotProfile`: wrappers for offline deepTools behavior.
- Modify `tests/data/e2e/fakebin/fake_bio_tool.py`: generate deterministic deepTools matrix/profile data and current AME text output.
- Modify `conf/test.config`: enable TSS QC in the offline end-to-end profile.
- Modify `README.md`: document report paths, fields, TSS definition, FRiP source, and missing-data behavior.
- Modify `CHANGELOG.md`: record the new consolidated dashboard.

### Task 1: QC data model, parsers, and calculations

**Files:**
- Create: `bin/qc_dashboard.py`
- Create: `tests/unit/test_qc_dashboard.py`

**Interfaces:**
- Consumes: UTF-8 JSON/TSV files already produced by demultiplexing, `LIBRARY_QC_CUSTOM`, `PEAK_QC`, `TSS_ENRICHMENT`, `MOTIF_QC_CUSTOM`, and AME.
- Produces:
  - `load_metadata(path: Path) -> dict[str, dict[str, object]]`
  - `read_demultiplex_metrics(paths: Sequence[Path]) -> dict[str, dict[str, object]]`
  - `read_library_metrics(paths: Sequence[Path]) -> dict[str, dict[str, object]]`
  - `read_peak_metrics(paths: Sequence[Path]) -> dict[str, dict[str, object]]`
  - `read_tss_profile(path: Path, *, before_bp: int = 3000, bin_size: int = 10) -> list[tuple[int, float]]`
  - `calculate_tss_enrichment(profile: Sequence[tuple[int, float]], *, flank_bp: int = 100) -> float | None`
  - `read_top_ame(path: Path, *, limit: int = 10) -> list[dict[str, object]]`
  - `build_report_data(metadata: Mapping[str, Mapping[str, object]], demultiplex: Mapping[str, Mapping[str, object]], libraries: Mapping[str, Mapping[str, object]], peaks: Mapping[str, Mapping[str, object]], tss: Mapping[str, Mapping[str, object]], motifs: Mapping[str, Mapping[str, object]], top_motifs: Mapping[str, Sequence[Mapping[str, object]]], annotation_status: str) -> dict[str, object]`
  - `DashboardInputError(ValueError)`

- [ ] **Step 1: Write failing metadata and table-parser tests**

Add tests which create temporary inputs and require strict identities:

```python
def test_load_metadata_rejects_duplicate_sample_ids(tmp_path):
    path = tmp_path / "metadata.json"
    path.write_text(json.dumps([
        {"sample_id": "S1", "library_id": "L1", "assay_target": "CTCF", "is_control": False},
        {"sample_id": "S1", "library_id": "L1", "assay_target": "CTCF", "is_control": False},
    ]))
    with pytest.raises(qc.DashboardInputError, match="duplicate sample_id S1"):
        qc.load_metadata(path)


def test_read_demultiplex_metrics_preserves_library_and_sample_counts(tmp_path):
    path = tmp_path / "demultiplex.metrics.json"
    path.write_text(json.dumps({
        "library_id": "L1",
        "total_reads": 100,
        "assigned_reads": 80,
        "ambiguous_reads": 5,
        "unassigned_reads": 15,
        "assignment_counts": {"S1": 30, "S2": 50},
    }))
    rows = qc.read_demultiplex_metrics([path])
    assert rows["L1"]["assigned_fraction"] == pytest.approx(0.8)
    assert rows["L1"]["assignment_counts"]["S1"] == 30
```

Also test blank identities, duplicate library IDs, duplicate metric names, negative
counts, non-finite numbers, and a peak-QC metric/value table containing `frip`.

- [ ] **Step 2: Run parser tests and verify the expected import failure**

Run:

```bash
python3 -m pytest tests/unit/test_qc_dashboard.py -k 'metadata or demultiplex or peak' -v
```

Expected: FAIL because `bin/qc_dashboard.py` does not exist.

- [ ] **Step 3: Implement strict reusable parsers**

Create an executable script with these constants and validation helpers:

```python
SCHEMA_VERSION = 1
TSS_BEFORE_BP = 3000
TSS_BIN_SIZE = 10
TSS_FLANK_BP = 100
TOP_MOTIF_LIMIT = 10


class DashboardInputError(ValueError):
    """Raised when a QC input cannot be joined without ambiguity."""


def finite_number(value: object, *, label: str, minimum: float | None = None) -> float:
    if isinstance(value, bool):
        raise DashboardInputError(f"{label} must be numeric")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise DashboardInputError(f"{label} must be numeric") from error
    if not math.isfinite(parsed) or (minimum is not None and parsed < minimum):
        raise DashboardInputError(f"{label} must be finite and >= {minimum}")
    return parsed
```

`load_metadata` must require `sample_id`, `library_id`, `assay_target`, and a
boolean `is_control`; preserve `input_group`, `control_id`, and
`expected_motif`; and return records sorted by sample ID.

`read_demultiplex_metrics` must calculate fractions from counts rather than
trust cached fractions. `read_library_metrics` reads one-row TSVs keyed by
`sample_id`. `read_peak_metrics` converts metric/value TSVs to one record and
requires a unique `sample_id`.

- [ ] **Step 4: Write failing TSS and AME tests**

Use a realistic deepTools profile table with label columns followed by 600
numeric bins. Require the parser to take the longest numeric suffix, assign
positions from `-3000` through `2990` in 10-bp increments, use index 300 as the center, and average
the first and last ten bins as the 200-bp combined flank:

```python
def test_tss_enrichment_uses_center_and_terminal_100bp(tmp_path):
    values = [2.0] * 600
    values[300] = 12.0
    path = tmp_path / "S1.tss_profile.tsv"
    path.write_text("sample\tgroup\t" + "\t".join(map(str, values)) + "\n")
    profile = qc.read_tss_profile(path)
    assert profile[0] == (-3000, 2.0)
    assert profile[300] == (0, 12.0)
    assert qc.calculate_tss_enrichment(profile) == pytest.approx(6.0)


def test_tss_zero_flank_returns_missing():
    profile = [(index * 10 - 3000, 0.0) for index in range(600)]
    assert qc.calculate_tss_enrichment(profile) is None
```

Add tests for the wrong bin count, `nan`, missing profiles, and top-ten AME
ordering. AME parsing must reuse `motif_qc.read_ame`, ignore the
`__NO_PEAKS__` sentinel, and order by adjusted significance then rank and motif
ID.

- [ ] **Step 5: Implement TSS and AME parsing**

`read_tss_profile` must ignore blank/comment-only lines and identify a row with
exactly `(2 * before_bp) // bin_size` trailing finite numeric fields. More than
one candidate data row is an incompatible schema and must fail.

`calculate_tss_enrichment` must:

```python
center_value = profile[len(profile) // 2][1]
flank_bins = TSS_FLANK_BP // TSS_BIN_SIZE
flank_values = [value for _, value in profile[:flank_bins]]
flank_values += [value for _, value in profile[-flank_bins:]]
flank_mean = statistics.fmean(flank_values)
return None if flank_mean == 0 else center_value / flank_mean
```

Import `read_ame` from the adjacent `motif_qc.py` and serialize each
`AmeRecord` as `motif_id`, `motif_alt_id`, `adjusted_p_value`, `p_value`,
`effect`, `positive_sequences`, and `rank`.

- [ ] **Step 6: Write and implement joined report-data tests**

Test one IgG plus one CTCF sample sharing a library. Require:

```python
assert data["samples"][0]["sample_id"] == "S_CTCF"
assert data["samples_by_id"]["S_CTCF"]["demultiplex"]["sample_assigned_reads"] == 50
assert data["samples_by_id"]["S_CTCF"]["peak"]["frip"] == pytest.approx(0.42)
assert data["samples_by_id"]["S_IgG"]["peak"]["frip"] is None
assert data["samples_by_id"]["S_IgG"]["motif"]["status"] == "not_applicable_control"
assert any("S_IgG" in warning["message"] for warning in data["warnings"])
```

`build_report_data` must initialize every optional family before overlaying
available data, join demultiplexing through `library_id`, distinguish
`target`/`control`, and reject metric records for unknown sample IDs. Missing
peak and motif data for controls is normal and should not receive a biological
failure label.

- [ ] **Step 7: Run unit tests and commit**

Run:

```bash
python3 -m pytest tests/unit/test_qc_dashboard.py -v
python3 -m pytest tests/unit/test_motif_qc.py -v
```

Expected: all tests pass.

Commit:

```bash
git add bin/qc_dashboard.py tests/unit/test_qc_dashboard.py
git commit -m "feat: aggregate consolidated QC metrics"
```

### Task 2: Deterministic TSV, JSON, SVG, and HTML report outputs

**Files:**
- Modify: `bin/qc_dashboard.py`
- Modify: `tests/unit/test_qc_dashboard.py`

**Interfaces:**
- Consumes: the `dict[str, object]` returned by `build_report_data`.
- Produces:
  - `write_qc_summary_tsv(data: Mapping[str, object], path: Path) -> None`
  - `write_qc_summary_json(data: Mapping[str, object], path: Path) -> None`
  - `write_top_motifs_tsv(data: Mapping[str, object], path: Path) -> None`
  - `write_tss_profiles_tsv(data: Mapping[str, object], path: Path) -> None`
  - `render_dashboard(data: Mapping[str, object]) -> str`
  - CLI output files `qc_dashboard.html`, `qc_summary.tsv`,
    `qc_summary.json`, `top_motifs.tsv`, and `tss_profiles.tsv`.

- [ ] **Step 1: Write failing serialization tests**

Define and assert this exact stable `QC_SUMMARY_COLUMNS` list containing
metadata, demultiplexing, alignment/library, peak/FRiP, TSS, and
expected-motif fields:

```python
QC_SUMMARY_COLUMNS = [
    "sample_id",
    "library_id",
    "input_group",
    "assay_target",
    "is_control",
    "control_id",
    "expected_motif",
    "total_read_pairs",
    "assigned_read_pairs",
    "ambiguous_read_pairs",
    "unassigned_read_pairs",
    "assigned_fraction",
    "ambiguous_fraction",
    "unassigned_fraction",
    "sample_assigned_reads",
    "sample_assignment_fraction",
    "raw_total_reads",
    "mapped_percent",
    "properly_paired_percent",
    "mapq_filtered_reads",
    "mapq_filtered_fragments",
    "mapq_filtered_fraction",
    "markdup_examined_reads",
    "duplicate_total",
    "duplicate_percent",
    "mitochondrial_percent",
    "estimated_library_size",
    "insert_size_total_pairs",
    "insert_size_min",
    "insert_size_q25",
    "insert_size_mean",
    "insert_size_median",
    "insert_size_q75",
    "insert_size_max",
    "peak_count",
    "total_covered_bases",
    "total_fragments",
    "fragments_in_peaks",
    "frip",
    "tss_status",
    "tss_enrichment",
    "expected_motif_status",
    "best_motif_id",
    "best_adjusted_p_value",
    "ame_status",
    "warning_count",
]
```

Require `None` to serialize as an empty TSV field and JSON `null`.

```python
def test_serializers_use_stable_columns_and_json_null(tmp_path, report_data):
    qc.write_outputs(report_data, tmp_path)
    rows = list(csv.DictReader(
        (tmp_path / "qc_summary.tsv").open(newline=""), delimiter="\t"
    ))
    assert list(rows[0]) == qc.QC_SUMMARY_COLUMNS
    assert rows[0]["frip"] == ""
    payload = json.loads((tmp_path / "qc_summary.json").read_text())
    assert payload["schema_version"] == 1
    assert payload["samples"][0]["peak"]["frip"] is None
```

Require `top_motifs.tsv` to contain no IgG rows and at most ten rows per target.
Require `tss_profiles.tsv` columns to be
`sample_id`, `position_bp`, `signal`.

- [ ] **Step 2: Implement deterministic serializers**

Write UTF-8 files with `newline=""`, tab delimiters, sample-ID ordering, motif
rank ordering, and TSS position ordering. JSON must use:

```python
json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
handle.write("\n")
```

Store metric definitions in JSON, including:

```python
"tss_enrichment": {
    "formula": "center_bin_signal / mean(terminal_100bp_flanks)",
    "before_bp": 3000,
    "after_bp": 3000,
    "bin_size_bp": 10,
    "center_position_bp": 0,
    "flank_bp_per_side": 100,
}
```

- [ ] **Step 3: Write failing self-contained HTML tests**

Require semantic section IDs and absence of remote resources:

```python
html = qc.render_dashboard(report_data)
for section_id in (
    "run-overview", "demultiplexing", "alignment",
    "peaks-frip", "tss-enrichment", "motif-enrichment",
):
    assert f'id="{section_id}"' in html
assert "<svg" in html
assert "https://" not in html
assert "http://" not in html
assert "overall pass" not in html.lower()
assert "overall fail" not in html.lower()
assert 'class="qc-pass"' not in html
assert 'class="qc-fail"' not in html
```

Also assert that warnings, `NA`, target/control legends, expected motifs, and
top AME rows are HTML-escaped and visible.

- [ ] **Step 4: Implement inline chart and report rendering**

Create small focused render helpers named `render_bar_chart`,
`render_line_chart`, `render_table`, and `render_dashboard`. The table helper
must use this escaping pattern, and the chart helpers must apply the same rule
to titles, labels, sample IDs, and tooltip text:

```python
def render_table(columns, rows, *, empty_message):
    if not rows:
        return f'<p class="empty">{html.escape(empty_message)}</p>'
    header = "".join(
        f"<th>{html.escape(label)}</th>" for _, label in columns
    )
    body = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(format_value(row.get(key)))}</td>"
            for key, _ in columns
        )
        + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"
```

Use `html.escape` on every external string. Use a fixed color palette with
solid target marks and outlined/hatched controls; never encode pass/fail.
Render responsive inline SVGs with view boxes, titles, axis labels, tooltips
using SVG `<title>`, and a text/table fallback under each chart. Include no
JavaScript dependency and no external URL.

- [ ] **Step 5: Implement the command-line entry point and atomic output**

Expose:

```text
qc_dashboard.py
  --metadata sample_metadata.json
  --demux-dir demux_inputs
  --library-dir library_inputs
  --insert-dir insert_inputs
  --peak-dir peak_inputs
  --tss-dir tss_inputs
  --motif-dir motif_inputs
  --ame-dir ame_inputs
  --annotation-status {skipped_no_annotation,computed_bed,computed_gtf}
  --outdir qc_dashboard
```

Create all five outputs in a temporary sibling directory and use
`os.replace()` only after every file validates as non-empty. On
`DashboardInputError`, print `qc_dashboard.py: error: <message>` and exit 2.

- [ ] **Step 6: Run the full dashboard unit suite and commit**

Run:

```bash
python3 -m pytest tests/unit/test_qc_dashboard.py -v
python3 -m compileall -q bin
```

Expected: all tests pass and compilation is silent.

Commit:

```bash
git add bin/qc_dashboard.py tests/unit/test_qc_dashboard.py
git commit -m "feat: render self-contained QC dashboard"
```

### Task 3: Export reproducible TSS profile data

**Files:**
- Modify: `modules/local/tss_enrichment.nf`
- Create: `tests/data/e2e/tss.bed`
- Create: `tests/data/e2e/fakebin/computeMatrix`
- Create: `tests/data/e2e/fakebin/plotProfile`
- Modify: `tests/data/e2e/fakebin/fake_bio_tool.py`
- Modify: `conf/test.config`
- Modify: `tests/integration/test_qc.sh`

**Interfaces:**
- Consumes: the existing RPKM bigWig and GTF/BED annotation.
- Produces: the existing TSS files plus
  `<sample_id>.tss_profile.tsv` from deepTools `plotProfile --outFileNameData`.
- Changes `TSS_ENRICHMENT.out.profiles` to:
  `tuple(meta, bed, matrix_gz, matrix_tsv, profile_png, profile_tsv, status_tsv)`.

- [ ] **Step 1: Add failing structural and fixture tests**

Update `tests/integration/test_qc.sh` to require:

```python
assert '--outFileNameData "${outputStem}.tss_profile.tsv"' in tss
assert 'path("*.tss_profile.tsv")' in tss
assert "meta, bed, matrix, matrixTable, profile, profileTable, status" in qc
```

Add `tests/data/e2e/tss.bed`:

```text
chrMini	199	200	tss_1	0	+
chrMini	799	800	tss_2	0	-
```

Set `params.tss_bed` to that file in `conf/test.config`.

- [ ] **Step 2: Run the QC integration test and verify failure**

Run:

```bash
bash tests/integration/test_qc.sh
```

Expected: FAIL because the TSS profile-data output and updated tuple do not yet
exist.

- [ ] **Step 3: Extend the TSS process**

Add:

```bash
plotProfile \
    --matrixFile "${outputStem}.tss_matrix.gz" \
    --outFileName "${outputStem}.tss_profile.png" \
    --outFileNameData "${outputStem}.tss_profile.tsv" \
    --plotTitle "TSS enrichment"
```

Add the profile TSV between the PNG and status in the output tuple and update
every destructuring expression in `subworkflows/local/qc.nf`.

- [ ] **Step 4: Add deterministic fake deepTools behavior**

The `computeMatrix` fake must write the requested matrix and matrix TSV. The
`plotProfile` fake must write a PNG placeholder plus a data row with 600 bins:
flanks equal `2`, center bin equal `12`, and all other bins equal `4`. Add both
tool names to `version()` and `handlers`.

Each wrapper must contain:

```bash
#!/usr/bin/env bash
exec "$(dirname "$0")/fake_bio_tool.py" --fake-tool "$(basename "$0")" "$@"
```

Also update `fake_ame` to emit AME `--text` TSV to stdout, matching the current
production module rather than requiring `--oc`.

- [ ] **Step 5: Run QC and environment checks, then commit**

Run:

```bash
bash tests/integration/test_qc.sh
bash tests/integration/check_envs.sh
```

Expected: both pass.

Commit:

```bash
git add modules/local/tss_enrichment.nf subworkflows/local/qc.nf \
  conf/test.config tests/data/e2e/tss.bed tests/data/e2e/fakebin/computeMatrix \
  tests/data/e2e/fakebin/plotProfile tests/data/e2e/fakebin/fake_bio_tool.py \
  tests/integration/test_qc.sh
git commit -m "feat: export TSS profile data"
```

### Task 4: Nextflow dashboard process and channel wiring

**Files:**
- Create: `modules/local/qc_dashboard.nf`
- Modify: `subworkflows/local/motifs.nf`
- Modify: `subworkflows/local/qc.nf`
- Modify: `main.nf`
- Modify: `tests/integration/test_qc.sh`
- Modify: `tests/integration/test_e2e.sh`

**Interfaces:**
- Consumes:
  - base64-encoded JSON list of validated sample metadata;
  - raw demultiplex JSON files;
  - library and insert-size custom TSVs;
  - peak metric TSVs and width histograms;
  - TSS average-profile and status TSVs;
  - expected-motif custom TSVs;
  - AME result directories and status TSVs;
  - one validated annotation-status value.
- Produces `QC_DASHBOARD.out` channels:
  - `report`: `qc_dashboard.html`
  - `summary_tsv`: `qc_summary.tsv`
  - `summary_json`: `qc_summary.json`
  - `top_motifs`: `top_motifs.tsv`
  - `tss_profiles`: `tss_profiles.tsv`
  - `versions`: `qc_dashboard_versions.yml`

- [ ] **Step 1: Add failing wiring assertions**

Update `tests/integration/test_qc.sh` to require:

```python
"include { QC_DASHBOARD } from '../../modules/local/qc_dashboard'" in qc
"QC_DASHBOARD(" in qc
"qc_dashboard_report = QC_DASHBOARD.out.report" in qc
"qc_summary_tsv = QC_DASHBOARD.out.summary_tsv" in qc
"known_motifs = AME.out.results" in motifs
"known_motif_statuses = AME.out.status" in motifs
```

Update `tests/integration/test_e2e.sh` so the static QC call expects eleven
arguments: the existing nine plus `motif_ame_results_ch` and
`motif_ame_status_ch`.

- [ ] **Step 2: Run structural tests and verify failure**

Run:

```bash
bash tests/integration/test_qc.sh
bash tests/integration/test_e2e.sh
```

Expected: structural assertions fail before any runtime assertions.

- [ ] **Step 3: Implement `QC_DASHBOARD`**

Use `label 'process_standard'`, `envs/python.yml`, and
`python:3.12.3-slim-bookworm`. Publish only the five report files to:

```groovy
publishDir "${params.outdir}/reports/qc_dashboard",
    mode: 'copy',
    overwrite: true,
    pattern: '{qc_dashboard.html,qc_summary.tsv,qc_summary.json,top_motifs.tsv,tss_profiles.tsv}'
```

Decode metadata with Python’s `base64` module, ensure every staged input
directory exists with `mkdir -p`, invoke `qc_dashboard.py` with the interface
from Task 2, and write:

```yaml
QC_DASHBOARD:
  python: Python 3.12.3
  qc_dashboard.py: repository
```

- [ ] **Step 4: Wire sample metadata and detailed artifacts in `QC`**

Add `motif_ame_results` and `motif_ame_statuses` to `take:`. Build metadata
from `safe_library_metrics`, sort by `sample_id`, reject duplicates, serialize
with `groovy.json.JsonOutput`, and base64-encode it before passing it as a
`val`.

Collect each file family into closed list values using the established
`reduce([files: []])` pattern. Collect both the TSS profile-data and status
files, peak histogram and metric files, raw demultiplex JSON, AME directories,
and AME status files. Invoke `QC_DASHBOARD` after the existing custom-format
processes and include its version record in `versions_ch`.

Emit all five dashboard outputs without changing the existing MultiQC emits.

- [ ] **Step 5: Expose AME channels and top-level outputs**

In `main.nf`, initialize:

```groovy
motif_ame_results_ch = Channel.empty()
motif_ame_status_ch = Channel.empty()
```

When motifs are enabled, assign `MOTIFS.out.known_motifs` and
`MOTIFS.out.known_motif_statuses`. Pass both to `QC`. Emit the dashboard HTML,
summary TSV/JSON, top motifs, and TSS profiles from `NANOCUT`.

Extend `WRITE_COMPLETION_SUMMARY` inputs with dashboard HTML and
`qc_summary.tsv`, validate both are non-empty, and add:

```text
qc_dashboard	reports/qc_dashboard/qc_dashboard.html
qc_summary	reports/qc_dashboard/qc_summary.tsv
```

- [ ] **Step 6: Run structural checks and Nextflow lint**

Run:

```bash
bash tests/integration/test_qc.sh
nextflow lint main.nf
```

Expected: QC assertions pass and lint reports no errors.

- [ ] **Step 7: Commit the workflow integration**

```bash
git add modules/local/qc_dashboard.nf subworkflows/local/motifs.nf \
  subworkflows/local/qc.nf main.nf tests/integration/test_qc.sh \
  tests/integration/test_e2e.sh
git commit -m "feat: integrate consolidated QC dashboard"
```

### Task 5: End-to-end verification and user documentation

**Files:**
- Modify: `tests/integration/test_e2e.sh`
- Modify: `tests/integration/test_docs.sh`
- Modify: `README.md`
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: the complete test-profile workflow.
- Produces: documented, verified report paths and machine-readable schemas.

- [ ] **Step 1: Add failing end-to-end dashboard assertions**

Require these files to be non-empty:

```python
results / "reports/qc_dashboard/qc_dashboard.html"
results / "reports/qc_dashboard/qc_summary.tsv"
results / "reports/qc_dashboard/qc_summary.json"
results / "reports/qc_dashboard/top_motifs.tsv"
results / "reports/qc_dashboard/tss_profiles.tsv"
```

Parse and assert:

```python
assert {row["sample_id"] for row in summary_rows} == {"MINI_IgG", "MINI_CTCF"}
ctcf = next(row for row in summary_rows if row["sample_id"] == "MINI_CTCF")
igg = next(row for row in summary_rows if row["sample_id"] == "MINI_IgG")
assert float(ctcf["frip"]) > 0
assert float(ctcf["tss_enrichment"]) == 6.0
assert igg["is_control"] == "true"
assert igg["frip"] == ""
assert igg["expected_motif_status"] == "not_applicable_control"
assert json_payload["schema_version"] == 1
assert len(top_motif_rows) == 1
assert top_motif_rows[0]["motif_alt_id"] == "CTCF"
assert len(tss_rows) == 1200
assert "https://" not in dashboard_html
```

Require the resumed run to cache `QC_DASHBOARD` alongside existing processes.

- [ ] **Step 2: Run end-to-end test and fix only fixture-contract defects**

Run:

```bash
bash tests/integration/test_e2e.sh
```

Expected: the first run publishes all dashboard outputs; the second run is
fully resumed/cached. If the runtime is unavailable, the script must still pass
its direct and structural checks and print its existing explicit skip message.

- [ ] **Step 3: Document reports and metric interpretation**

Update the README output tree and QC sections to document:

- MultiQC remains at `reports/multiqc/multiqc_report.html`;
- the detailed report is `reports/qc_dashboard/qc_dashboard.html`;
- `qc_summary.tsv`, `qc_summary.json`, `top_motifs.tsv`, and
  `tss_profiles.tsv` are reusable data products;
- the original FRiP value remains at
  `qc/peaks/<sample_id>/<sample_id>.peak_qc.tsv`, row `frip`;
- TSS enrichment uses position 0 divided by the mean first/last 100 bp;
- IgG controls are visually separated and have no expected-motif result;
- missing optional outputs appear as `NA` warnings rather than zero;
- the dashboard is descriptive and applies no biological thresholds.

Add those phrases to `tests/integration/test_docs.sh`.

- [ ] **Step 4: Update the changelog**

Under the current unreleased/added section, record the self-contained
dashboard, joined TSV/JSON summaries, TSS scalar/profile export, and top-ten AME
summary without claiming biological classification.

- [ ] **Step 5: Run the complete repository verification**

Run:

```bash
python3 -m pytest tests/unit -v
bash tests/integration/check_envs.sh
bash tests/integration/test_docs.sh
bash tests/integration/test_demultiplex.sh
bash tests/integration/test_alignment.sh
bash tests/integration/test_peaks.sh
bash tests/integration/test_motifs.sh
bash tests/integration/test_qc.sh
bash tests/integration/test_e2e.sh
nextflow lint main.nf
git diff --check
```

Expected: all available tests pass, Nextflow lint reports no errors, and
`git diff --check` is silent.

- [ ] **Step 6: Commit documentation and final integration assertions**

```bash
git add README.md CHANGELOG.md tests/integration/test_docs.sh \
  tests/integration/test_e2e.sh
git commit -m "docs: describe consolidated QC reports"
```

### Task 6: Final review and pull-request readiness

**Files:**
- Review all files changed since `origin/agent/bulk-nanocut-tag-pipeline`.

**Interfaces:**
- Consumes: completed Tasks 1–5.
- Produces: a verified feature branch ready to push and open as a pull request.

- [ ] **Step 1: Review scope and accidental changes**

Run:

```bash
git status --short
git diff --stat origin/agent/bulk-nanocut-tag-pipeline...HEAD
git diff --name-status origin/agent/bulk-nanocut-tag-pipeline...HEAD
```

Expected: only the design, plan, dashboard implementation, workflow wiring,
fixtures/tests, README, and changelog are present. `.superpowers/` remains
untracked and is not staged.

- [ ] **Step 2: Re-run the completion verification**

Use the `superpowers:verification-before-completion` skill and re-run the full
Task 5 verification from a clean shell. Record the exact pass/skip output for
the PR description.

- [ ] **Step 3: Review generated HTML with the miniature fixture**

Open the generated `qc_dashboard.html` and verify that:

- all six sections render;
- target and IgG styling are distinguishable;
- tables remain readable at narrow width;
- missing values show `NA`;
- no chart labels overlap;
- the file opens with networking disabled.

- [ ] **Step 4: Prepare the pull request**

Summarize:

- the hybrid MultiQC plus detailed-dashboard architecture;
- new report/data paths;
- the TSS score definition;
- target/control handling;
- missing-data behavior;
- all verification commands and any runtime skips.

Do not push or create the pull request until the user requests publication or
the active execution workflow reaches its explicit publication handoff.
