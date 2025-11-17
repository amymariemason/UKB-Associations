#### Utility helpers for working with UKB outcome control sheet outputs
#### (originally based on code from Scott Richie)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, tidyverse) 

# Simple infix helper to provide a default when a value is NULL.  Keeps
# the curate_settings return value backwards compatible with the naming
# used elsewhere in the repository.
`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

# Function to load a dataset from the RAP into a local directory.  A
# sensible default work_dir is provided so the function can also be used
# locally without repeatedly specifying the target folder.
load_from_rap <- function(rap_path, work_dir = "./Inputs") {
  if (!dir.exists(sprintf("%s/input_data/", work_dir))) {
    system(sprintf("mkdir -p %s/input_data", work_dir), ignore.stdout = FALSE)
  }
  fname <- basename(rap_path)
  if (!file.exists(sprintf("%s/input_data/%s", work_dir, fname))) {
    system(sprintf("dx download '%s' -o '%s/input_data/%s'", rap_path, work_dir, fname), ignore.stdout = FALSE)
  }
  fread(sprintf("%s/input_data/%s", work_dir, fname), na.strings = c("", "NA"))
}

# test: Load ICD codes associated with cancer register
# cancer_icd10 <- load_from_rap("common/Cancer Register/icd10_codes.csv", "./Inputs")

# Helper to determine whether a code vector contains any usable entries.
# This mirrors the behaviour of the settings objects produced by
# `parse_control_sheet()` in Read_settings.R: code vectors are either
# empty character vectors or vectors of regex patterns.
has_codes <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(FALSE)
  }
  any(!is.na(x) & trimws(as.character(x)) != "")
}

# Summarise the definition list produced by parse_control_sheet() into a
# tidy table indicating whether each outcome has codes for each coding
# system (ICD9, ICD10, Read codes, etc.).  This makes it easy to review
# newly generated settings objects before running the full pipeline.
summarise_outcome_definitions <- function(outcome_defs) {
  if (is.null(outcome_defs) || length(outcome_defs) == 0) {
    stop("No outcome definitions supplied for summarising.")
  }

  summaries <- lapply(outcome_defs, function(def) {
    data.table(
      outcome_id = def$outcome_id,
      label = def$label,
      prevalent_incident = isTRUE(def$prevalent_incident),
      use_primary_care = isTRUE(def$use_primary_care),
      has_icd9 = has_codes(def$icd9),
      has_icd10 = has_codes(def$icd10),
      has_death = has_codes(def$death),
      has_self_report = any(c(
        has_codes(def$self_report_20002),
        has_codes(def$self_report_20004),
        has_codes(def$self_report_6150),
        has_codes(def$self_report_6152),
        has_codes(def$self_report_6153),
        has_codes(def$self_report_6177)
      )),
      has_procedures = any(c(has_codes(def$procedures_opcs3), has_codes(def$procedures_opcs4))),
      has_cancer_registry = has_codes(def$cancer_registry),
      has_primary_care_read = any(c(has_codes(def$read_v2), has_codes(def$read_v3))),
      match_icd10_primary_care = isTRUE(def$match_icd10_to_read)
    )
  })

  rbindlist(summaries, fill = TRUE)
}

# Produce a single-row overview for each code type indicating whether any
# outcome definition uses it, alongside the list of definitions requiring
# that code type.  Intended to sit alongside the row-wise outcome summary
# from `summarise_outcome_definitions()` for quick review.
data_requirements <- function(definitions_summary) {
  if (!is.data.table(definitions_summary)) {
    stop("definitions_summary must be a data.table from summarise_outcome_definitions().")
  }

  code_flags <- list(
    prevalent_incident="prevalent_incident",
    icd9 = "has_icd9",
    icd10 = "has_icd10",
    death = "has_death",
    self_report = "has_self_report",
    procedures = "has_procedures",
    cancer_registry = "has_cancer_registry",
    primary_care_read = "has_primary_care_read",
    match_icd10_primary_care = "match_icd10_primary_care"
  )

  # 1. Create the named tibble of any_required = TRUE/FALSE
  any_required_tbl <- tibble(
    !!! setNames(
      map(code_flags, ~ any(definitions_summary[[.x]])),
      names(code_flags)
    )
  )
  
  # 2. Create a named list of outcomes per code type
  outcomes_list <- map(
    code_flags,
    ~ definitions_summary[get(.x) == TRUE, outcome_id]
  )
  
  # Name the list for clarity
  names(outcomes_list) <- names(code_flags)
  
  # Return both objects together
  list(
    datasets = any_required_tbl,
    outcome_lists = outcomes_list
  )
}

# Given a settings object produced by parse_control_sheet(), return a
# curated version with an attached summary table.  
curate_settings <- function(settings) {
  if (is.null(settings$outcomes_def)) {
    stop("The settings object does not contain an 'outcomes_def' element.")
  }

  summary_tbl <- summarise_outcome_definitions(settings$outcomes_def)
  data_requirements <- data_requirements(summary_tbl)

  list(
    filename_out = settings$filename_out,
    matched_outcomes = settings$matched_outcomes,
    diabetes_flag = settings$diabetes_flag,
    outcomes_def = settings$outcomes_def,
    data_requirements = data_requirements$datasets,
    data_loop_lists = data_requirements$outcome_lists
  )
}
