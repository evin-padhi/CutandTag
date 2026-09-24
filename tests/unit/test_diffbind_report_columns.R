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
  chrom = c("chr1", "chr2"),
  start = c(10L, 20L),
  end = c(30L, 45L),
  interval_key = c("chr1:11:30", "chr2:21:45"),
  stringsAsFactors = FALSE
)
one_based_report <- data.frame(
  Chr = c("chr1", "chr2"), Start = c(11L, 21L), End = c(30L, 45L),
  stringsAsFactors = FALSE
)
normalized_one_based <- normalize_diffbind_report_coordinates(
  one_based_report, bed_intervals, "test report"
)
stopifnot(
  identical(normalized_one_based$coordinate_mode, "1-based report coordinates"),
  identical(normalized_one_based$report$interval_key, bed_intervals$interval_key)
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

cat("PASS: DiffBind report fields and coordinate joins are validated.\n")
