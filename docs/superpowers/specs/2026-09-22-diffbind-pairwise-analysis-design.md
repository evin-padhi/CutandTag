# DiffBind Pairwise Analysis Design

## Objective

Add a downstream differential binding analysis to the bulk nano-CUT&Tag Nextflow pipeline. For every assay with enough target samples and conditions, use DiffBind to test all pairwise condition contrasts over the existing condition-aware, blacklist-filtered consensus narrow-peak intervals. Run two analyses over the same intervals: target-only counts and matched-IgG-subtracted counts. Save combined peak-level results, a compact comparison summary table, and clean diagnostic plots for each contrast.

## Current pipeline context

The pipeline already:

- Reads one manifest row per sample, with `condition`, `assay_target`, `is_control`, and `control_id` fields.
- Aligns paired-end reads and produces filtered BAMs.
- Calls matched-IgG narrow peaks for every target sample.
- Merges reproducible narrow peaks within each condition, then merges these condition sets into a consensus BED. An optional blacklist is removed from this BED.
- Counts paired-end fragments in the consensus intervals and writes a raw fragment count matrix.

The new analysis starts after the consensus BED is made. It does not change the existing peak calls or fragment count matrix.

## Scientific model and assumptions

- The assay under test is identified by `assay_target`; control samples are never conditions in the differential model.
- `condition` identifies the biological groups. Each target sample is one replicate. The initial supported model is an unpaired `~ Condition` design because the current manifest does not define matched blocks or repeated measures.
- Each condition must have at least two target samples for inferential testing. If an assay has fewer than two conditions, or a condition has fewer than two replicates, the differential analysis for that assay must stop with an actionable message. The pipeline must not silently treat technical libraries as biological replicates.
- Each target's `control_id` resolves to its matched IgG sample. The same target and IgG BAMs are used in both counting modes.
- The pipeline has no spike-in normalization input. Results describe relative enrichment differences across the measured libraries. They do not establish absolute genome-wide occupancy changes.
- The planned dataset has three conditions (`disomy`, `ts18`, and `ts21`), so it will produce three pairwise contrasts for each assay.

## Design decision

Use DiffBind with the explicit consensus BED as its fixed peak set. Disable any automatic peak calling and automatic blacklist/greylist filtering inside DiffBind. This keeps all samples and both counting modes on the same genomic intervals and avoids hidden changes to the tested regions.

Run `dba.count` and the DESeq2 analysis twice per assay, with identical sample metadata, BAMs, consensus intervals, fragment settings, and normalization settings:

1. **Target only:** Count target BAM fragments in the consensus intervals. Do not attach or subtract controls in this mode.
2. **Target minus matched IgG:** Attach each sample's matched IgG BAM and apply DiffBind's scaled control subtraction. Keep the same intervals and DESeq2 settings as target-only mode.

Make control scaling and subtraction explicit in the implementation. Do not rely on DiffBind defaults, because defaults can depend on whether a blacklist or greylist is present. Confirm the exact DiffBind API options and their behavior against the pinned package version during implementation.

Use DESeq2 for both modes. For the three-condition dataset, test all unordered pairs: `disomy_vs_ts18`, `disomy_vs_ts21`, and `ts18_vs_ts21`. Label every result with its assay, contrast, count mode, peak ID, and genomic coordinates. Calculate Benjamini-Hochberg adjusted p-values separately within each contrast and count mode. Use adjusted p-value < 0.05 as the default significance threshold and include that threshold and adjustment scope in output metadata.

## Data flow

1. Group target samples by assay and condition from the normalized manifest. Resolve each target's `control_id` to the IgG BAM and metadata.
2. Collect the assay's filtered target BAMs, required matched IgG BAMs, and its consensus BED.
3. Build a DiffBind sample sheet in the task from typed inputs. Preserve sample IDs and condition labels. Do not place task-local file paths in a workflow-scope generated file.
4. Count fragments over the fixed consensus BED in target-only and target-minus-IgG modes.
5. Run the three pairwise DESeq2 contrasts in each mode.
6. From normalized target-only counts over the consensus peaks, calculate replicate correlations for each assay and condition.
7. Combine the differential results and create the summary table and plots.

If any target lacks its declared IgG BAM, if BAMs cannot be mapped unambiguously to manifest sample IDs, or if the consensus BED has no valid intervals, fail with a clear message before statistical testing.

## Outputs

Write outputs under `differential_binding/<assay_target>/`:

- `diffbind_results.tsv`: one row per tested peak, contrast, and count mode. Include coordinates, both condition means, log2 fold change, raw p-value, adjusted p-value, count mode, and significance flag. Preserve the DiffBind/DESeq2 values needed to reproduce summary totals.
- `diffbind_comparison_summary.tsv`: one row per contrast and count mode, with condition names, sample counts, tested peak count, significant peak count, counts with positive and negative log2 fold change, and FDR threshold.
- One MA plot and one volcano plot per contrast. Each plot presents target-only and target-minus-IgG results as clearly labelled facets or panels. Do not add plot titles or subtitles. Use consistent axes and significance styling across contrasts.
- `replicate_correlation.tsv`: one row per target-sample pair within the same assay and condition, with sample IDs, condition, number of consensus peaks used, Pearson correlation, and Spearman correlation. Calculate correlations from normalized target-only counts after a `log2(count + 1)` transformation. Do not calculate correlations from IgG-subtracted counts.
- One correlation heatmap per assay using all target samples, with samples annotated by condition. Add scatter plots for target-sample pairs within each condition. Plot transformed normalized target-only peak counts, use equal axis scales within each scatter plot, and label only points that have an extreme contribution to the correlation.
- A run log and software-version record for the DiffBind R package and its R/Bioconductor runtime.

The current raw fragment count matrix remains available and unchanged. These new outputs provide the statistical analysis; the raw matrix remains useful for auditing counts.

## Error handling and reporting

- Validate the experimental design before launching DiffBind: at least two conditions and at least two target replicates per condition are required for every tested assay.
- Require each declared matched control BAM to be present and readable for IgG-subtracted mode.
- Stop with the offending assay, sample, condition, or contrast in the error message when an input is invalid.
- If a valid contrast yields no peaks below the FDR threshold, write the results and summary row with zero significant peaks. Do not treat this as a pipeline error.
- Report both count modes separately. Do not merge their p-values or imply that they are interchangeable estimators. Adjust p-values within each contrast and mode, as specified above.
- Treat replicate correlations as descriptive QC. Do not use them to remove samples automatically or to claim that replicates are biologically valid. If correlation cannot be calculated because values are constant or too few values are finite, report it as unavailable with a reason.

## Scope exclusions

- No new peak caller or peak-union strategy. The existing consensus BED is the fixed testing universe.
- No paired, blocked, batch-adjusted, or covariate model until the manifest defines those design variables and the experiment supports them.
- No automatic selection of a preferred count mode. Users compare both modes and interpret IgG-subtracted estimates with the raw-target analysis as context.
- No claim of absolute occupancy change without a suitable external normalization strategy such as spike-in calibration.

## Acceptance conditions

1. Every eligible assay uses the same consensus intervals and sample set in both count modes.
2. Each target sample is assigned to exactly one manifest condition and its declared matched IgG is used only for the subtraction mode.
3. All pairwise contrasts are produced for each eligible assay, with labels that state the reference and comparison conditions clearly.
4. Results TSV contains peak coordinates, effect estimates, p-values, adjusted p-values, count mode, and contrast labels.
5. Summary TSV reports per-contrast and per-mode test and significance counts.
6. Every contrast has an MA plot and volcano plot that distinguish count modes.
7. Replicate correlations use normalized, transformed target-only counts over the same consensus intervals used for differential testing, and include both Pearson and Spearman measures.
8. Correlation heatmaps and within-condition scatter plots identify sample conditions and do not add titles or subtitles.
9. Invalid replication, control mapping, or interval inputs fail before statistical testing with actionable messages.
10. The README documents the method, inputs, output files, interpretation, replicate QC, and relative-binding limitation.

## References

- DiffBind package page: <https://bioconductor.org/packages/release/bioc/html/DiffBind.html>
- DiffBind user manual: <https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf>
- DiffBind vignette: <https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf>
