#!/usr/bin/env Rscript

script_path <- file.path("bin", "diffbind_analysis.R")
script_expressions <- parse(script_path)
function_definition <- Filter(
  function(expression) {
    is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      identical(expression[[2]], as.name("resolve_concentration_columns"))
  },
  as.list(script_expressions)
)
if (length(function_definition) != 1) {
  stop("Expected one resolve_concentration_columns() definition in the analysis script.")
}
test_environment <- new.env(parent = baseenv())
eval(function_definition[[1]], envir = test_environment)
resolve_concentration_columns <- test_environment$resolve_concentration_columns

coordinate_definition <- Filter(
  function(expression) {
    is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      identical(expression[[2]], as.name("normalize_diffbind_report_coordinates"))
  },
  as.list(script_expressions)
)
if (length(coordinate_definition) != 1) {
  stop("Expected one normalize_diffbind_report_coordinates() definition.")
}
eval(coordinate_definition[[1]], envir = test_environment)
normalize_diffbind_report_coordinates <- test_environment$normalize_diffbind_report_coordinates

binding_counts_definition <- Filter(
  function(expression) {
    is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      identical(expression[[2]], as.name("extract_full_sample_counts"))
  },
  as.list(script_expressions)
)
if (length(binding_counts_definition) != 1) {
  stop("Expected one extract_full_sample_counts() definition.")
}
eval(binding_counts_definition[[1]], envir = test_environment)
extract_full_sample_counts <- test_environment$extract_full_sample_counts

legacy_report <- data.frame(
  Conc = 10,
  Conc_group1 = 8,
  Conc_group2 = 12,
  check.names = FALSE
)
stopifnot(identical(
  resolve_concentration_columns(legacy_report, "disomy", "trisomy18"),
  c("Conc_group1", "Conc_group2")
))

condition_named_report <- data.frame(
  Conc = 10,
  Conc_disomy = 8,
  Conc_trisomy18 = 12,
  check.names = FALSE
)
stopifnot(identical(
  resolve_concentration_columns(condition_named_report, "disomy", "trisomy18"),
  c("Conc_disomy", "Conc_trisomy18")
))

unexpected_report <- data.frame(Chr = "chr1", Fold = 1, check.names = FALSE)
error_message <- tryCatch(
  resolve_concentration_columns(unexpected_report, "disomy", "trisomy18"),
  error = conditionMessage
)
stopifnot(
  grepl("Could not identify DiffBind group concentration columns", error_message),
  grepl("Chr, Fold", error_message)
)

bed_intervals <- data.frame(
  chrom = c("chr1", "chr2", "chr3", "chr4"),
  start = c(10L, 20L, 30L, 40L),
  end = c(30L, 45L, 55L, 70L),
  interval_key = c("chr1:11:30", "chr2:21:45", "chr3:31:55", "chr4:41:70"),
  stringsAsFactors = FALSE
)
one_based_report <- data.frame(
  Chr = c("chr1", "chr2", "chr3", "chr4"),
  Start = c(11L, 21L, 31L, 41L),
  End = c(30L, 45L, 55L, 70L),
  stringsAsFactors = FALSE
)
normalized_one_based <- normalize_diffbind_report_coordinates(
  one_based_report, bed_intervals, "test report"
)
stopifnot(
  identical(normalized_one_based$coordinate_mode, "1-based report coordinates"),
  identical(normalized_one_based$report$interval_key, bed_intervals$interval_key)
)

partial_report <- data.frame(
  Chr = c("chr1", "chr3", "chr4"),
  Start = c(11L, 31L, 41L),
  End = c(30L, 55L, 70L),
  Fold = c(1.2, -0.8, 0.4),
  stringsAsFactors = FALSE
)
normalized_partial <- normalize_diffbind_report_coordinates(
  partial_report, bed_intervals, "partial test report"
)
partial_join <- merge(
  bed_intervals,
  normalized_partial$report[, c("interval_key", "Fold")],
  by = "interval_key",
  all.x = TRUE,
  sort = FALSE
)
stopifnot(
  nrow(normalized_partial$report) == 3,
  nrow(partial_join) == nrow(bed_intervals),
  sum(is.na(partial_join$Fold)) == 1
)

binding_matrix <- data.frame(
  CHR = rev(bed_intervals$chrom),
  START = rev(bed_intervals$start + 1L),
  END = rev(bed_intervals$end),
  sample_a = c(4, 3, 2, 1),
  sample_b = c(40, 30, 20, 10),
  stringsAsFactors = FALSE
)
full_sample_counts <- extract_full_sample_counts(
  binding_matrix, bed_intervals, c("sample_a", "sample_b")
)
stopifnot(
  nrow(full_sample_counts) == nrow(bed_intervals),
  setequal(full_sample_counts$interval_key, bed_intervals$interval_key),
  full_sample_counts$sample_a[match("chr1:11:30", full_sample_counts$interval_key)] == 1,
  full_sample_counts$sample_b[match("chr4:41:70", full_sample_counts$interval_key)] == 40
)

zero_based_report <- transform(one_based_report, Start = Start - 1L)
normalized_zero_based <- normalize_diffbind_report_coordinates(
  zero_based_report, bed_intervals, "test report"
)
stopifnot(
  identical(
    normalized_zero_based$coordinate_mode,
    "0-based report coordinates normalized to 1-based BED keys"
  ),
  identical(normalized_zero_based$report$interval_key, bed_intervals$interval_key)
)

unmatched_report <- transform(one_based_report, Chr = "other")
coordinate_error <- tryCatch(
  normalize_diffbind_report_coordinates(unmatched_report, bed_intervals, "test report"),
  error = conditionMessage
)
stopifnot(
  grepl("coordinates do not match the consensus intervals", coordinate_error),
  grepl("direct matches=0", coordinate_error)
)

duplicate_report <- one_based_report[c(1, 2, 3, 1), ]
duplicate_error <- tryCatch(
  normalize_diffbind_report_coordinates(duplicate_report, bed_intervals, "duplicate test report"),
  error = conditionMessage
)
stopifnot(grepl("duplicate interval coordinates", duplicate_error))

cat("PASS: DiffBind report fields and coordinate joins are validated.\n")
