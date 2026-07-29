# Motif BED Column Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Allow motif background preparation to combine MACS2 peak files and BED3 blacklists without BEDTools field-count failures.

**Architecture:** Convert each exclusion source to chromosome, start, and end coordinates before concatenation. Keep the conversion inside `PREPARE_MOTIF_SEQUENCES` so input files and downstream interfaces remain unchanged.

**Tech Stack:** Nextflow DSL2, Bash, AWK, BEDTools, shell integration tests

## Global Constraints

- Preserve the original staged peak and blacklist inputs.
- Preserve all existing motif outputs and channel interfaces.
- Exclusion sorting and merging consume only valid BED3 rows.

---

### Task 1: Normalize motif background exclusions

**Files:**
- Modify: `tests/integration/test_motifs.sh`
- Modify: `modules/local/prepare_motif_sequences.nf`

**Interfaces:**
- Consumes: staged `final.peaks.bed` and optional `blacklist/regions.bed`
- Produces: `background_exclusions.unsorted.bed` containing consistent BED3 rows

- [ ] **Step 1: Write the failing regression assertion**

Add assertions requiring AWK BED3 projection of both final peaks and the
optional blacklist before either source reaches
`background_exclusions.unsorted.bed`.

- [ ] **Step 2: Run the focused motif test and verify failure**

Run: `bash tests/integration/test_motifs.sh`

Expected: FAIL because the current implementation copies nine-column final
peaks and appends a three-column blacklist without normalization.

- [ ] **Step 3: Implement BED3 normalization**

Use AWK to print `chromosome`, `start`, and `end` from each non-comment,
nonblank input row. Write normalized final peaks to the exclusion file and
append normalized blacklist coordinates only when a blacklist is supplied.

- [ ] **Step 4: Run focused and repository checks**

Run:

```bash
bash tests/integration/test_motifs.sh
for test_script in tests/integration/test_*.sh; do
    bash "$test_script"
done
nextflow lint .
```

Expected: all available tests pass; unavailable optional runtimes report their
existing explicit skip messages.

- [ ] **Step 5: Commit and publish**

Stage only the design, plan, test, and motif module. Commit the regression fix,
push `codex/fix-motif-bed-columns`, and open a pull request against the remote
default branch.
