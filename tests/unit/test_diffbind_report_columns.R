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

cat("PASS: DiffBind report concentration columns support legacy and condition-based names.\n")
