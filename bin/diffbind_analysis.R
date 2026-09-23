#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DiffBind)
  library(DESeq2)
  library(tidyverse)
})

log_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

parse_args <- function(args) {
  if (length(args) %% 2 != 0) {
    stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
  }
  keys <- args[seq(1, length(args), by = 2)]
  values <- args[seq(2, length(args), by = 2)]
  if (any(!str_detect(keys, "^--[a-z][a-z-]*$")) || anyDuplicated(keys)) {
    stop("Argument names must be unique values such as --sample-sheet.", call. = FALSE)
  }
  set_names(as.list(values), str_remove(keys, "^--"))
}

required_arg <- function(options, name) {
  value <- options[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop("Required argument is missing: --", str_replace_all(name, "_", "-"), call. = FALSE)
  }
  value
}

check_bam <- function(path, sample_id, role) {
  if (is.na(path) || !nzchar(path) || !file.exists(path) || file.access(path, 4) != 0) {
    stop(role, " BAM for sample ", sample_id, " is missing or unreadable: ", path, call. = FALSE)
  }
  index_candidates <- c(paste0(path, ".bai"), str_replace(path, "\\.bam$", ".bai"))
  if (!any(file.exists(index_candidates))) {
    stop(role, " BAM index for sample ", sample_id, " is missing beside ", path, call. = FALSE)
  }
}

validate_bed <- function(path) {
  bed <- suppressWarnings(read_tsv(
    path,
    col_names = FALSE,
    comment = "#",
    show_col_types = FALSE,
    progress = FALSE,
    col_types = cols(.default = col_character())
  ))
  if (nrow(bed) == 0 || ncol(bed) < 3) {
    stop("Consensus BED must contain at least one interval and three columns.", call. = FALSE)
  }
  starts <- suppressWarnings(as.numeric(bed[[2]]))
  ends <- suppressWarnings(as.numeric(bed[[3]]))
  if (anyNA(bed[[1]]) || any(!nzchar(bed[[1]])) || anyNA(starts) || anyNA(ends) ||
      any(!is.finite(starts)) || any(!is.finite(ends)) ||
      any(starts != floor(starts)) || any(ends != floor(ends)) ||
      any(starts > .Machine$integer.max) || any(ends > .Machine$integer.max) ||
      any(starts < 0) || any(ends <= starts)) {
    stop("Consensus BED has invalid chromosome, start, or end values.", call. = FALSE)
  }
  tibble(
    chrom = bed[[1]],
    start = as.integer(starts),
    end = as.integer(ends)
  ) %>%
    transmute(
      peak_id = sprintf("peak_%07d", row_number()),
      chrom,
      start,
      end,
      interval_key = paste(chrom, start + 1L, end, sep = ":")
    )
}

validate_samples <- function(path) {
  samples <- read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("", "NA")
  ) %>%
    mutate(across(c(sample_id, condition, bam_reads, control_id, bam_control), as.character))

  needed <- c("sample_id", "condition", "bam_reads", "control_id", "bam_control")
  missing <- setdiff(needed, names(samples))
  if (length(missing) > 0) {
    stop("Sample sheet is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (nrow(samples) == 0 || anyNA(samples$sample_id) || any(!nzchar(samples$sample_id))) {
    stop("Sample sheet must contain target sample IDs.", call. = FALSE)
  }
  if (anyDuplicated(samples$sample_id)) {
    stop("Sample sheet contains duplicate target sample IDs.", call. = FALSE)
  }
  if (anyNA(samples$condition) || any(!nzchar(samples$condition))) {
    stop("Every target sample must have a condition.", call. = FALSE)
  }
  if (any(!file.exists(samples$bam_reads))) {
    bad <- samples %>% filter(!file.exists(bam_reads)) %>% slice(1)
    stop("Target BAM is missing for sample ", bad$sample_id, ": ", bad$bam_reads, call. = FALSE)
  }
  walk2(samples$bam_reads, samples$sample_id, ~check_bam(.x, .y, "Target"))
  control_rows <- samples %>% filter(!is.na(control_id), !is.na(bam_control))
  if (nrow(control_rows) != nrow(samples)) {
    bad <- samples %>% filter(is.na(control_id) | is.na(bam_control)) %>% slice(1)
    stop("Matched IgG BAM mapping is missing for target sample ", bad$sample_id, call. = FALSE)
  }
  if (any(!file.exists(samples$bam_control))) {
    bad <- samples %>% filter(!file.exists(bam_control)) %>% slice(1)
    stop("Matched IgG BAM is missing for target sample ", bad$sample_id, ": ", bad$bam_control, call. = FALSE)
  }
  walk2(samples$bam_control, samples$sample_id, ~check_bam(.x, .y, "Matched IgG"))

  condition_counts <- samples %>% count(condition, name = "replicates")
  if (nrow(condition_counts) < 2) {
    stop("Differential binding needs at least two conditions.", call. = FALSE)
  }
  insufficient <- condition_counts %>% filter(replicates < 2)
  if (nrow(insufficient) > 0) {
    stop(
      "At least two target replicates are required in every condition. Insufficient: ",
      assay, ": ",
      paste0(insufficient$condition, " (", insufficient$replicates, ")", collapse = ", "),
      call. = FALSE
    )
  }
  samples %>% arrange(condition, sample_id)
}

make_dba <- function(samples, assay, count_mode, consensus_bed, sample_peak_file) {
  use_control <- identical(count_mode, "target_minus_igg")
  sheet <- samples %>%
    transmute(
      SampleID = sample_id,
      Tissue = assay,
      Factor = assay,
      Condition = condition,
      Replicate = as.integer(as.factor(sample_id)),
      bamReads = bam_reads,
      Peaks = sample_peak_file,
      PeakCaller = "bed"
    )
  if (use_control) {
    sheet <- sheet %>%
      mutate(
        bamControl = samples$bam_control,
        ControlID = samples$control_id
      )
  }

  log_message("Starting DiffBind counts for ", assay, " in mode ", count_mode)
  dba_object <- dba(
    sampleSheet = as.data.frame(sheet),
    config = list(
      AnalysisMethod = DBA_DESEQ2,
      doBlacklist = FALSE,
      doGreylist = FALSE,
      RunParallel = FALSE,
      bCorPlot = FALSE,
      singleEnd = FALSE,
      fragments = FALSE,
      inter.feature = TRUE
    ),
    peakCaller = "bed",
    bRemoveM = FALSE,
    bRemoveRandom = FALSE
  )
  dba_object <- dba.count(
    dba_object,
    peaks = consensus_bed,
    summits = FALSE,
    filter = 0,
    minCount = 0,
    bRemoveDuplicates = FALSE,
    bScaleControl = use_control,
    bSubControl = use_control,
    bUseSummarizeOverlaps = TRUE,
    bParallel = FALSE
  )
  dba_object <- dba.contrast(
    dba_object,
    design = "~Condition",
    minMembers = 2,
    categories = DBA_CONDITION
  )
  dba_object <- dba.analyze(
    dba_object,
    method = DBA_DESEQ2,
    design = "~Condition",
    bBlacklist = FALSE,
    bGreylist = FALSE,
    bParallel = FALSE
  )
  dba_object
}

safe_contrast <- function(value) {
  str_replace_all(value, "[^A-Za-z0-9._-]", "_")
}

report_for_contrast <- function(dba_object, contrast_id, contrast_row, count_mode, bed, sample_ids, assay, fdr) {
  report <- dba.report(
    dba_object,
    contrast = contrast_id,
    method = DBA_DESEQ2,
    th = 1,
    bUsePval = FALSE,
    bCounts = FALSE,
    bNormalized = TRUE,
    DataType = DBA_DATA_FRAME
  ) %>%
    as.data.frame(check.names = FALSE) %>%
    as_tibble()
  count_report <- dba.report(
    dba_object,
    contrast = contrast_id,
    method = DBA_DESEQ2,
    th = 1,
    bUsePval = FALSE,
    bCounts = TRUE,
    bNormalized = FALSE,
    DataType = DBA_DATA_FRAME
  ) %>%
    as.data.frame(check.names = FALSE) %>%
    as_tibble()

  required_columns <- c("Chr", "Start", "End", "Conc_group1", "Conc_group2", "Fold", "p-value", "FDR")
  absent <- setdiff(required_columns, names(report))
  if (length(absent) > 0) {
    stop("DiffBind report is missing expected fields: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  count_columns <- intersect(sample_ids, names(count_report))
  if (length(count_columns) == 0) {
    stop("DiffBind report has no per-sample counts for contrast ", contrast_id, call. = FALSE)
  }
  report_stats <- report %>%
    transmute(
      interval_key = paste(Chr, Start, End, sep = ":"),
      mean_group1 = as.numeric(Conc_group1),
      mean_group2 = as.numeric(Conc_group2),
      log2_fold_change = as.numeric(Fold),
      p_value = as.numeric(`p-value`),
      adjusted_p_value = as.numeric(FDR)
    )
  result <- bed %>%
    left_join(report_stats, by = "interval_key") %>%
    mutate(
      assay_target = assay,
      contrast = paste0(contrast_row$Group1, "_vs_", contrast_row$Group2),
      condition1 = as.character(contrast_row$Group1),
      condition2 = as.character(contrast_row$Group2),
      count_mode = count_mode,
      fdr_threshold = fdr,
      significant = !is.na(adjusted_p_value) & adjusted_p_value < fdr
    )

  if (nrow(result) != nrow(bed)) {
    stop("Internal error while retaining consensus intervals in contrast ", contrast_id, call. = FALSE)
  }
  counts <- count_report %>%
    transmute(
      interval_key = paste(Chr, Start, End, sep = ":"),
      across(all_of(count_columns))
    )
  list(results = result, counts = counts)
}

write_diffbind_results <- function(samples, bed, assay, fdr, outdir, consensus_bed) {
  modes <- c("target_only", "target_minus_igg")
  all_results <- list()
  target_counts <- NULL
  summary_rows <- list()
  expected_pairs <- combn(sort(unique(samples$condition)), 2, simplify = FALSE)
  sample_peak_file <- file.path(outdir, "diffbind_sample_peaks.bed")
  bed %>%
    transmute(chrom, start, end, peak_id, score = 1L) %>%
    write_tsv(sample_peak_file, col_names = FALSE)

  for (count_mode in modes) {
    dba_object <- make_dba(samples, assay, count_mode, consensus_bed, sample_peak_file)
    contrasts <- as_tibble(dba.show(dba_object, bContrasts = TRUE))
    if (!all(c("Group1", "Group2") %in% names(contrasts))) {
      stop("DiffBind did not report condition pair names.", call. = FALSE)
    }
    selected <- contrasts %>%
      mutate(contrast_id = row_number()) %>%
      mutate(pair_key = map2_chr(Group1, Group2, ~paste(sort(c(.x, .y)), collapse = "\t"))) %>%
      filter(pair_key %in% map_chr(expected_pairs, ~paste(sort(.x), collapse = "\t")))
    if (nrow(selected) != length(expected_pairs)) {
      stop("DiffBind did not create every pairwise condition contrast.", call. = FALSE)
    }

    for (contrast_id in seq_len(nrow(selected))) {
      contrast_row <- selected[contrast_id, ]
      log_message("Analyzing ", assay, " ", count_mode, " contrast ", contrast_row$Group1, " vs ", contrast_row$Group2)
      extracted <- report_for_contrast(
        dba_object,
        contrast_row$contrast_id,
        contrast_row,
        count_mode,
        bed,
        samples$sample_id,
        assay,
        fdr
      )
      result <- extracted$results
      all_results[[length(all_results) + 1]] <- result %>%
        select(
          peak_id, chrom, start, end, assay_target,
          condition1, condition2, contrast, count_mode,
          mean_group1, mean_group2, log2_fold_change, p_value,
          adjusted_p_value, fdr_threshold, significant
        )
      if (count_mode == "target_only") {
        target_counts[[length(target_counts) + 1]] <- extracted$counts
      }
      summary_rows[[length(summary_rows) + 1]] <- result %>%
        summarise(
          condition1 = first(condition1),
          condition2 = first(condition2),
          contrast = first(contrast),
          count_mode = first(count_mode),
          sample_count1 = sum(samples$condition == first(condition1)),
          sample_count2 = sum(samples$condition == first(condition2)),
          tested_peaks = n(),
          significant_peaks = sum(significant),
          positive_log2_fold_change = sum(significant & log2_fold_change > 0, na.rm = TRUE),
          negative_log2_fold_change = sum(significant & log2_fold_change < 0, na.rm = TRUE),
          fdr_threshold = first(fdr_threshold)
        )

    }
  }

  results <- bind_rows(all_results) %>% arrange(count_mode, contrast, chrom, start, end)
  summaries <- bind_rows(summary_rows) %>% arrange(contrast, count_mode)
  write_tsv(results, file.path(outdir, "diffbind_results.tsv"), na = "")
  write_tsv(summaries, file.path(outdir, "diffbind_comparison_summary.tsv"), na = "")
  plot_data <- results %>%
    mutate(
      count_mode = recode(count_mode,
        target_only = "Target only",
        target_minus_igg = "Target minus matched IgG"
      ),
      minus_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
      significance = if_else(significant, "FDR < threshold", "Not significant")
    )
  for (contrast in unique(plot_data$contrast)) {
    contrast_data <- plot_data %>% filter(.data$contrast == contrast)
    prefix <- safe_contrast(contrast)
    ma_plot <- ggplot(contrast_data, aes(x = (mean_group1 + mean_group2) / 2, y = log2_fold_change, color = significance)) +
      geom_point(size = 0.7, alpha = 0.65) +
      geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
      facet_wrap(vars(count_mode), scales = "fixed") +
      scale_color_manual(values = c("FDR < threshold" = "#B23A48", "Not significant" = "#7A8793")) +
      labs(x = "Mean normalized signal (log2)", y = "Log2 fold change", color = NULL) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank())
    volcano_plot <- ggplot(contrast_data, aes(x = log2_fold_change, y = minus_log10_p, color = significance)) +
      geom_point(size = 0.7, alpha = 0.65) +
      geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
      geom_hline(yintercept = -log10(fdr), color = "grey55", linetype = "dashed", linewidth = 0.3) +
      facet_wrap(vars(count_mode), scales = "fixed") +
      scale_color_manual(values = c("FDR < threshold" = "#B23A48", "Not significant" = "#7A8793")) +
      labs(x = "Log2 fold change", y = "-log10(raw p-value)", color = NULL) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank())
    ggsave(file.path(outdir, paste0(prefix, "_ma.pdf")), ma_plot, width = 8, height = 4.6)
    ggsave(file.path(outdir, paste0(prefix, "_volcano.pdf")), volcano_plot, width = 8, height = 4.6)
  }
  list(results = results, target_counts = bind_rows(target_counts), summary = summaries)
}

write_replicate_qc <- function(samples, bed, target_counts, outdir) {
  raw_counts <- target_counts %>%
    pivot_longer(cols = all_of(samples$sample_id), names_to = "sample_id", values_to = "normalized_count") %>%
    filter(!is.na(normalized_count)) %>%
    distinct(interval_key, sample_id, .keep_all = TRUE)
  count_wide <- raw_counts %>%
    select(interval_key, sample_id, normalized_count) %>%
    pivot_wider(names_from = sample_id, values_from = normalized_count, values_fill = 0) %>%
    arrange(interval_key)
  raw_matrix <- count_wide %>% select(all_of(samples$sample_id)) %>% as.matrix()
  storage.mode(raw_matrix) <- "numeric"
  if (!any(raw_matrix > 0)) {
    stop("Target-only consensus-peak counts are all zero; replicate correlations cannot be computed.", call. = FALSE)
  }
  size_factors <- estimateSizeFactorsForMatrix(raw_matrix, type = "poscounts")
  normalized_matrix <- sweep(raw_matrix, 2, size_factors, "/")
  count_long <- as_tibble(normalized_matrix, .name_repair = "minimal") %>%
    set_names(samples$sample_id) %>%
    mutate(interval_key = count_wide$interval_key, .before = 1) %>%
    pivot_longer(
      cols = all_of(samples$sample_id),
      names_to = "sample_id",
      values_to = "normalized_count"
    ) %>%
    mutate(log2_count = log2(pmax(as.numeric(normalized_count), 0) + 1))
  peak_index <- bed %>% select(peak_id, interval_key)
  count_long <- count_long %>% left_join(peak_index, by = "interval_key")
  pairs <- samples %>%
    group_by(condition) %>%
    group_modify(~ {
      sample_pairs <- combn(.x$sample_id, 2, simplify = FALSE)
      tibble(
        sample_id_1 = map_chr(sample_pairs, 1),
        sample_id_2 = map_chr(sample_pairs, 2)
      )
    }) %>%
    ungroup()
  correlations <- pmap_dfr(pairs, function(condition, sample_id_1, sample_id_2) {
    pair_values <- count_long %>%
      filter(sample_id %in% c(sample_id_1, sample_id_2)) %>%
      select(peak_id, sample_id, log2_count) %>%
      pivot_wider(names_from = sample_id, values_from = log2_count) %>%
      filter(is.finite(.data[[sample_id_1]]), is.finite(.data[[sample_id_2]]))
    x <- pair_values[[sample_id_1]]
    y <- pair_values[[sample_id_2]]
    reason <- case_when(
      length(x) < 2 ~ "fewer than two finite peak values",
      n_distinct(x) < 2 || n_distinct(y) < 2 ~ "one sample has constant signal across peaks",
      TRUE ~ NA_character_
    )
    tibble(
      condition = condition,
      sample_id_1 = sample_id_1,
      sample_id_2 = sample_id_2,
      peaks_used = length(x),
      pearson = if (is.na(reason)) cor(x, y, method = "pearson") else NA_real_,
      spearman = if (is.na(reason)) cor(x, y, method = "spearman") else NA_real_,
      unavailable_reason = reason
    )
  })
  write_tsv(correlations, file.path(outdir, "replicate_correlation.tsv"), na = "")

  sample_matrix <- count_long %>%
    select(peak_id, sample_id, log2_count) %>%
    pivot_wider(names_from = sample_id, values_from = log2_count) %>%
    select(all_of(samples$sample_id)) %>%
    as.matrix()
  rownames(sample_matrix) <- count_long %>% distinct(peak_id) %>% pull(peak_id)
  correlation_matrix <- cor(sample_matrix, method = "pearson", use = "pairwise.complete.obs")
  correlation_long <- as_tibble(as.table(correlation_matrix), .name_repair = "minimal") %>%
    set_names(c("sample_id_1", "sample_id_2", "pearson")) %>%
    left_join(samples %>% select(sample_id, condition) %>% rename(condition_1 = condition), by = c("sample_id_1" = "sample_id")) %>%
    left_join(samples %>% select(sample_id, condition) %>% rename(condition_2 = condition), by = c("sample_id_2" = "sample_id")) %>%
    mutate(
      sample_id_1_label = paste0(sample_id_1, "\n", condition_1),
      sample_id_2_label = paste0(sample_id_2, "\n", condition_2)
    )
  heatmap <- ggplot(correlation_long, aes(sample_id_1_label, sample_id_2_label, fill = pearson)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(limits = c(-1, 1), low = "#3B6FB6", mid = "white", high = "#B23A48", midpoint = 0) +
    labs(x = NULL, y = NULL, fill = "Pearson r") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(outdir, "replicate_correlation_heatmap.pdf"), heatmap, width = 7, height = 6)

  for (row in seq_len(nrow(pairs))) {
    condition <- pairs$condition[[row]]
    sample_id_1 <- pairs$sample_id_1[[row]]
    sample_id_2 <- pairs$sample_id_2[[row]]
    scatter_data <- count_long %>%
      filter(sample_id %in% c(sample_id_1, sample_id_2)) %>%
      select(peak_id, sample_id, log2_count) %>%
      pivot_wider(names_from = sample_id, values_from = log2_count) %>%
      filter(is.finite(.data[[sample_id_1]]), is.finite(.data[[sample_id_2]])) %>%
      mutate(extreme = abs(as.numeric(scale(.data[[sample_id_1]] - .data[[sample_id_2]]))) > 3)
    scatter <- ggplot(scatter_data, aes(.data[[sample_id_1]], .data[[sample_id_2]])) +
      geom_abline(slope = 1, intercept = 0, color = "grey60", linewidth = 0.3) +
      geom_point(size = 0.7, alpha = 0.55, color = "#356A7A") +
      geom_text(data = filter(scatter_data, extreme), aes(label = peak_id), size = 2, check_overlap = TRUE, vjust = -0.5) +
      coord_equal() +
      labs(x = sample_id_1, y = sample_id_2) +
      theme_minimal(base_size = 9) +
      theme(panel.grid.minor = element_blank())
    ggsave(
      file.path(outdir, paste0("replicate_scatter_", safe_contrast(condition), "_", safe_contrast(sample_id_1), "_vs_", safe_contrast(sample_id_2), ".pdf")),
      scatter,
      width = 5,
      height = 5
    )
  }
  correlations
}

main <- function() {
  cli_options <- parse_args(commandArgs(trailingOnly = TRUE))
  sample_sheet <- required_arg(cli_options, "sample-sheet")
  consensus_bed <- required_arg(cli_options, "consensus-bed")
  assay <- required_arg(cli_options, "assay")
  outdir <- required_arg(cli_options, "outdir")
  fdr <- as.numeric(cli_options[["fdr"]] %||% "0.05")
  if (length(fdr) != 1 || is.na(fdr) || fdr <= 0 || fdr >= 1) {
    stop("--fdr must be greater than 0 and less than 1.", call. = FALSE)
  }
  if (!file.exists(sample_sheet) || file.access(sample_sheet, 4) != 0) {
    stop("Sample sheet is missing or unreadable: ", sample_sheet, call. = FALSE)
  }
  if (!file.exists(consensus_bed) || file.access(consensus_bed, 4) != 0) {
    stop("Consensus BED is missing or unreadable: ", consensus_bed, call. = FALSE)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(outdir, "diffbind_analysis.log")
  sink(log_file, split = TRUE)
  on.exit(sink(), add = TRUE)
  sink(log_file, type = "message", append = TRUE)
  on.exit(sink(type = "message"), add = TRUE)
  log_message("Validating inputs for assay ", assay)
  samples <- validate_samples(sample_sheet)
  bed <- validate_bed(consensus_bed)
  write_tsv(bed %>% select(-interval_key), file.path(outdir, "consensus_peak_ids.tsv"))
  log_message("Validated ", nrow(samples), " targets across ", n_distinct(samples$condition), " conditions and ", nrow(bed), " intervals")
  result <- write_diffbind_results(samples, bed, assay, fdr, outdir, consensus_bed)
  write_replicate_qc(samples, bed, result$target_counts, outdir)
  write_lines(c(
    "DIFFBIND_ANALYSIS:",
    paste0("  R: ", getRversion()),
    paste0("  DiffBind: ", packageVersion("DiffBind")),
    paste0("  DESeq2: ", packageVersion("DESeq2")),
    paste0("  tidyverse: ", packageVersion("tidyverse"))
  ), file.path(outdir, "diffbind_analysis_versions.yml"))
  log_message("Wrote DiffBind, comparison summary, replicate correlation tables, and plots to ", outdir)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
main()
