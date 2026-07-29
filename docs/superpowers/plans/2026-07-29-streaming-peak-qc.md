# Streaming Peak QC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prevent large BEDPE inputs from exhausting PEAK_QC memory by processing fragments incrementally.

**Architecture:** Introduce `iter_fragments(path)` for line-oriented parsing, make fragment normalization and name deduplication lazy, and count overlaps during the same traversal. Retain the current peak index and output interfaces.

**Tech Stack:** Python 3.12, unittest, Nextflow DSL2

## Global Constraints

- Repeated fragment names count once.
- FRiP, per-peak counts, filenames, and output schemas remain unchanged.
- Peak inputs may remain in memory because they are orders of magnitude smaller than BEDPE inputs.
- The CLI must not materialize the BEDPE fragment stream.

---

### Task 1: Stream BEDPE fragments through peak QC

**Files:**
- Modify: `tests/unit/test_peak_qc.py`
- Modify: `bin/peak_qc.py`

**Interfaces:**
- Consumes: `Iterable[Fragment | Sequence[object]]` and BEDPE paths
- Produces: unchanged JSON, metrics TSV, width histogram, and fragments-per-peak TSV outputs

- [ ] **Step 1: Write the failing laziness regression**

Add a guarded generator that yields one fragment, then asserts overlap
processing happened before it yields the second fragment. Invoke
`write_qc_outputs` with this generator and assert both fragments are counted.

- [ ] **Step 2: Verify the regression fails**

Run:

```bash
python3 -m unittest \
  tests.unit.test_peak_qc.PeakQcUnitTests.test_write_qc_outputs_processes_fragments_incrementally
```

Expected: failure because `write_qc_outputs` currently consumes the complete
generator while constructing normalized fragment lists.

- [ ] **Step 3: Implement lazy parsing and overlap counting**

Create `iter_fragments(path)` that opens the BEDPE file and yields one validated
`Fragment` at a time. Convert normalization and name deduplication helpers to
yield iterators, update overlap counting to increment its own total, and have
the CLI pass `iter_fragments(arguments.fragments)` directly to
`write_qc_outputs`.

- [ ] **Step 4: Verify focused and complete suites**

Run:

```bash
python3 -m unittest tests.unit.test_peak_qc
for test_script in tests/integration/test_*.sh; do
    bash "$test_script"
done
python3 -m unittest discover -s tests/unit -p 'test_*.py'
```

Expected: all available tests pass; optional Nextflow runtime checks may emit
their existing explicit skip messages.

- [ ] **Step 5: Commit and publish**

Stage the two implementation files and the streaming design/plan documents.
Commit, push `codex/fix-motif-bed-columns`, and create a pull request against
the repository default branch.
