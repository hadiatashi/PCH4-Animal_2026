# Quick software check for the simulation workflow.
# Runs the main script with B_AGC=50 in a separate output directory.
# This is not intended as the full bootstrap analysis.

main <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) == 1L) {
    runner_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    runner_dir <- normalizePath(getwd())
  }

  src <- file.path(runner_dir, "ch4_AGC_ME_indices_simulation_test.R")
  if (!file.exists(src)) {
    stop("Cannot find main script: ", src)
  }

  quick_dir <- file.path(runner_dir, "quick_test_output")
  dir.create(quick_dir, showWarnings = FALSE, recursive = TRUE)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(quick_dir)
  Sys.setenv(B_AGC = "50")

  cat("Quick-test directory: ", quick_dir, "\n", sep = "")
  run_env <- new.env(parent = globalenv())
  sys.source(src, envir = run_env, chdir = FALSE)
}

main()
