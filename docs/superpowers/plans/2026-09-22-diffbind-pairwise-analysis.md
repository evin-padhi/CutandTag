# DiffBind Pairwise Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add pairwise DiffBind analysis and descriptive replicate correlation QC for target samples in the bulk nano-CUT&Tag pipeline.

**Architecture:** Keep consensus peak creation unchanged. Extend the differential-binding subworkflow to map each target BAM to its manifest condition and matched IgG BAM, then run one isolated R task per assay on the shared consensus BED. The R program will run target-only and scaled IgG-subtracted DiffBind count modes, test all pairwise conditions with DESeq2, calculate target-only replicate correlations, and create the specified TSV and plot outputs.

**Tech Stack:** Nextflow DSL2, R, DiffBind, DESeq2, tidyverse, ggplot2, Bioconductor.

**Spec:** `docs/superpowers/specs/2026-09-22-diffbind-pairwise-analysis-design.md`

## Global Constraints

- Use the same blacklist-filtered consensus intervals for every sample, contrast, and count mode.
- Use matched IgG controls only for the IgG-subtracted count mode; use target-only normalized counts for replicate correlations.
- Adjust p-values with Benjamini-Hochberg separately within each contrast and count mode; use adjusted p-value < 0.05 as the default threshold.
- Use a `~ Condition` unpaired model. Require at least two target replicates per condition for inferential testing.
- Treat replicate correlations as descriptive QC. Do not remove samples automatically.
- Do not claim absolute genome-wide occupancy change without spike-in normalization.
- Use tidyverse syntax in R scripts. Plots must not have titles or subtitles and must use a clean, minimal style.
- Record task progress in log messages and capture R, DiffBind, and relevant Bioconductor package versions.
- Do not run local Docker builds for smoke tests. If an existing GitHub Actions workflow is suitable, use it for the pipeline smoke check after pushing.
- Do not change or remove the existing raw fragment-count matrix.

## Review Focus

- A target points to a missing, non-IgG, or wrong-assay `control_id`; cover this in manifest-to-BAM pairing validation and ensure subtraction stops before DiffBind.
- A target BAM or peak file is missing, duplicated, or joined to the wrong sample ID; validate exact one-to-one mappings before writing the task sample sheet.
- A condition has one target replicate or fewer; validate design eligibility before testing and report assay, condition, and sample count.
- A consensus BED is empty, malformed, or has invalid coordinates; reject it before counting and report the assay.
- Correlations encounter constant or non-finite vectors; write unavailable values with a reason and keep other sample-pair results.

---

## File Structure

- Create `envs/diffbind.yml`: pin an R/Bioconductor runtime and DiffBind version with DESeq2, tidyverse, and plotting dependencies.
- Create `bin/diffbind_analysis.R`: validate task-local sample and interval inputs; run both count modes and all contrasts; calculate correlation QC; write TSV, plots, logs, and package versions.
- Create `modules/local/diffbind_analysis.nf`: stage one assay's target BAMs, matched IgG BAMs, fixed consensus BED, and metadata; create the R sample sheet in task scope; run the R program and publish outputs.
- Modify `subworkflows/local/differential_binding.nf`: retain target metadata and controls, resolve `control_id` to IgG BAMs, pass complete assay inputs to DiffBind, and emit its results/logs/versions.
- Modify `main.nf`: expose DiffBind result outputs and include the new version records in version collection.
- Modify `README.md`: document eligibility rules, fixed consensus regions, the two count modes, outputs, correlation calculation, and limits on interpretation.
- Modify `CHANGELOG.md`: record pairwise differential binding and replicate QC.

## Task 1: Implement the DiffBind analysis program and pinned environment

**Files:**
- Create: `envs/diffbind.yml`
- Create: `bin/diffbind_analysis.R`

**Interfaces:**
- Consumes a task-local DiffBind sample sheet with target sample IDs, conditions, target BAM paths, and matched control BAM paths; a fixed BED file; and an output directory.
- Produces `diffbind_results.tsv`, `diffbind_comparison_summary.tsv`, `replicate_correlation.tsv`, per-contrast MA and volcano plots, per-assay correlation heatmap, within-condition replicate scatter plots, a log, and a software-version YAML file.

- [ ] Pin an R and Bioconductor release combination that supports the selected DiffBind version. Include only needed packages: DiffBind, DESeq2, tidyverse, and any direct dependencies required for plot rendering. Use strict `conda-forge` then `bioconda` channel order.
- [ ] Define named CLI arguments for sample sheet, consensus BED, assay label, FDR cutoff, and output directory. Reject missing files, unreadable BAMs, duplicate sample IDs, invalid condition labels, missing matched controls in subtraction mode, invalid BED rows, and empty intervals with actionable messages.
- [ ] Validate the design before count/test operations: at least two conditions and at least two target samples per condition. Use sample IDs in any validation message.
- [ ] Build one DiffBind object for target-only counts and one for scaled matched-IgG subtraction. Pass the consensus BED as the explicit fixed peak set and disable automatic peak calling and automatic blacklist/greylist filtering. Set subtraction and scaling options explicitly. Verify option names and behavior against the pinned DiffBind API before coding those calls.
- [ ] Use the same target records, intervals, fragment-count settings, and DESeq2 settings for both modes. Test every unordered pair of conditions. Label the contrast direction in all outputs and use BH adjustment within each contrast and mode.
- [ ] Extract per-peak coordinates, condition means, log2 fold change, raw p-value, adjusted p-value, count mode, contrast, and significance flag. Write one combined result table and one summary row per contrast and mode, including sample counts and significant positive/negative effect counts.
- [ ] Calculate Pearson and Spearman correlations from normalized target-only consensus-peak counts after `log2(count + 1)`. Write pairwise metrics for samples within each assay and condition. For constant or insufficient vectors, write unavailable values and a reason.
- [ ] Create clean ggplot2 MA and volcano plots per contrast with count modes shown in labelled panels; create a condition-annotated sample correlation heatmap per assay and scatter plots for replicate pairs within conditions. Use minimal themes and no titles/subtitles. Ensure all figures use the required stable names and readable sample labels.
- [ ] Log each validation, count, contrast, QC, and output step. Write package and R versions to a YAML file.

## Task 2: Add the isolated Nextflow task and wire target/control BAM inputs

**Files:**
- Create: `modules/local/diffbind_analysis.nf`
- Modify: `subworkflows/local/differential_binding.nf`

**Interfaces:**
- Consumes the assay-grouped target BAMs, their condition and control ID metadata, matched IgG BAMs, and the assay consensus BED from existing workflow outputs.
- Produces typed assay-level result, QC, log, and version channels for the top-level workflow.

- [ ] Define a process with the pinned `envs/diffbind.yml` environment and an appropriate versioned container if repository conventions require both Conda and container declarations.
- [ ] Accept metadata and BAM files as explicit task inputs. Build the sample sheet inside the task after Nextflow has localized the BAM paths; safely quote or serialize sample values without placing task-local paths in a workflow-scope file.
- [ ] In the subworkflow, retain the normalized metadata required to link each target's `control_id` to exactly one IgG BAM. Do not treat IgG controls as conditions or include them in target sample correlations.
- [ ] Validate that all target BAMs, peak calls, conditions, and control mappings have exact sample-ID matches before the analysis process runs.
- [ ] Join the assay consensus BED to the assay input group without changing the existing consensus merge or raw featureCounts matrix path.
- [ ] Emit all DiffBind result tables and plot files, plus log and versions. Use a stable output directory `differential_binding/<assay_target>/` and avoid collisions with existing outputs.

## Task 3: Expose outputs and document the analysis

**Files:**
- Modify: `main.nf`
- Modify: `README.md`
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes the new DiffBind subworkflow channels.
- Produces top-level workflow outputs for the result TSVs and figures, includes DiffBind versions in version collection, and documents the user-facing files and interpretation.

- [ ] Pass the existing target/control metadata, filtered BAM channels, and consensus BED into the updated differential-binding subworkflow.
- [ ] Add emitted DiffBind results, plots, logs, and version records to the top-level workflow output and version collection using existing naming patterns.
- [ ] Document all output names and locations, the two count modes, all pairwise comparisons, minimum replication rule, correlation inputs and statistics, and the absence of spike-in calibration.
- [ ] Update the changelog with the new analysis outputs.
- [ ] Review the final diff against every acceptance condition in the design spec and ensure no local Docker build or test run is added to the plan.

## Execution Notes

- Keep changes within the current PR worktree and branch.
- Follow existing Nextflow module and publication patterns in `modules/local/count_differential_fragments.nf` and `subworkflows/local/differential_binding.nf`.
- The sample sheet is an internal task input generated after BAM localization. It is not a new user-supplied manifest.
- Run no tests unless the user asks to test or verify the implementation. Do not claim the feature is validated by a full pipeline run unless such a run is explicitly requested and completed.
