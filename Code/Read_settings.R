
# Control sheet parsing utilities for the UK Biobank outcome pipeline
#
# These helpers are responsible for reading the bespoke outcome control
# sheet used by the legacy Stata workflow and transforming it into a
# structured list of outcome definitions that can be consumed by the R
# implementation.  The code expects the spreadsheet layout that ships
# with this repository:
#   * a `Settings` sheet describing run time configuration
#   * a `Definitions` sheet containing one row per outcome
#   * optional Read code mapping sheets (`Read_V2`, `Read_V3`)
#

library(readxl)
library(tidyverse)

# Normalise column names from the Excel workbook into lower snake case so
# we can address them reliably regardless of minor formatting
# differences.
normalise_colname <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  tolower(trimws(x))
}

# Convert "YES"/"NO" (case insensitive) or TRUE/FALSE like values into a
# logical flag.  Anything that cannot be interpreted defaults to FALSE.
parse_yes_no <- function(x) {
  if (is.logical(x)) {
    return(ifelse(is.na(x), FALSE, x))
  }
  x <- toupper(trimws(as.character(x)))
  ifelse(x %in% c("Y", "YES", "TRUE", "1"), TRUE, FALSE)
}

# Transform a comma separated string of codes into a vector of regular
# expressions.  The control sheet uses the convention where a trailing
# ".X" should behave as a wildcard prefix match.  Periods elsewhere in
# the code are stripped in line with the Stata implementation.
parse_code_vector <- function(x, wildcard = "[A-Z0-9]*") {
  if (is.null(x) || length(x) == 0) {
    return(character())
  }
  x <- x[!is.na(x)]
  x <- trimws(x)
  x <- x[x != "" & x != "-"]
  if (length(x) == 0) {
    return(character())
  }
  split_values <- unlist(strsplit(x, "[,]"))
  split_values <- trimws(split_values)
  split_values <- split_values[split_values != ""]
  if (length(split_values) == 0) {
    return(character())
  }
  split_values <- toupper(split_values)
  split_values <- gsub("\\.X$", wildcard, split_values, perl = TRUE)
  split_values <- gsub("\\.", "", split_values, fixed = FALSE)
  paste0("^", split_values, "$")
}

# Helper to read the Settings sheet.  The sheet stores configuration in a
# simple key/value layout and a list of outcomes to produce.

#test:  path<- "Inputs/bespoke_outcome_v3.xls"
read_settings_sheet <- function(path, CVD=F, CANCER=F) {
  settings_raw <- readxl::read_excel(path, sheet = "Settings", .name_repair = "minimal")
  if (nrow(settings_raw) == 0) {
    stop("The Settings sheet is empty – please check the control sheet.")
    }
  colnames(settings_raw) <- normalise_colname(colnames(settings_raw))
  if (!"setting" %in% colnames(settings_raw) || !"value" %in% colnames(settings_raw)) {
    stop("The Settings sheet must contain 'Setting' and 'Value' columns.")
    }
  
  #
  settings_raw <- settings_raw %>% filter(!is.na(setting))
  
  #get filename
  filename_out<- settings_raw %>% filter(setting=="Filename") %>% pull(value)
  
  if (CVD){
    filename_out<-"CVD"
    }
  if (CANCER){
    filename_out<-"Cancer"
    }
  
  filename_out<-ifelse(is.na(filename_out)|filename_out=="" ,"temp", filename_out)
  filename_out<-paste0(filename_out,"_",format(Sys.Date(), "%M%Y"))
                       
  # get diabetes
  diabetes_flag<- parse_yes_no(settings_raw %>% filter(setting=="diabetes") %>% pull(value))
                       
   # get list of outcomes wanted
   value_map <- settings_raw$value
   if (CVD){
    value_map<-settings_raw$cvd}
   if (CANCER){
    value_map<-settings_raw$cancer}
    
  value_map<-settings_raw$value
  value_map<- vapply(value_map, parse_yes_no, FUN.VALUE = logical(1))
  
  requested_outcomes <- character()
  
  keep_rows <- which(value_map==TRUE)
  requested_outcomes <- settings_raw$setting[keep_rows]
  requested_outcomes <- requested_outcomes[requested_outcomes != ""]
  
  return(list(
    filename= filename_out,
    outcomes = requested_outcomes,
    diabetes_flag = diabetes_flag
      )
    )
}

# Parse the Definitions sheet into a tidy data frame with one row per
# outcome and list columns holding the relevant code sets.
read_outcome_definitions <- function(path, outcome_filter = NULL) {
    defs_raw <- readxl::read_excel(path, sheet = "Definitions", .name_repair = "minimal")
    if (nrow(defs_raw) == 0) {
      stop("The Definitions sheet is empty – please check the control sheet.")
        }
   colnames(defs_raw) <- normalise_colname(colnames(defs_raw))
   defs_raw <- defs_raw[!is.na(defs_raw$suggested_variable_name) &
       trimws(defs_raw$suggested_variable_name) != "", ]
   if (!is.null(outcome_filter) && length(outcome_filter) > 0) {
   match_ids <- tolower(trimws(defs_raw$suggested_variable_name)) %in%
       tolower(trimws(outcome_filter))
       defs_raw <- defs_raw[match_ids, ]
     }
   if (nrow(defs_raw) == 0) {
       stop("No matching outcomes found in Definitions sheet for the requested set.")
   }
#1. Define the function to create a single list (row) of outcome parameters
    make_def <- function(row) {
      list(
        outcome_id = trimws(as.character(row[["suggested_variable_name"]])),
        label = trimws(as.character(row[["outcome_name"]])),
        prevalent_incident = parse_yes_no(row[["prevalent_incident"]]),
        use_primary_care = parse_yes_no(row[["use_primary_care"]]),
        match_icd10_to_read = parse_yes_no(row[["match_icd10_primary_care_codes"]]),
        icd9 = parse_code_vector(row[["icd_9_codes"]], wildcard = "[0-9]*"),
        icd10 = parse_code_vector(row[["icd_10_codes"]], wildcard = "[0-9]*"),
        death = parse_code_vector(row[["death_40001_40002"]], wildcard = "[0-9]*"),
        self_report_20002 = parse_code_vector(row[["self_report_20002"]], wildcard = "[0-9]*"),
        self_report_20004 = parse_code_vector(row[["self_report_20004"]], wildcard = "[0-9]*"),
        self_report_6150 = parse_code_vector(row[["self_report_6150"]], wildcard = "[0-9]*"),
        self_report_6152 = parse_code_vector(row[["self_report_6152"]], wildcard = "[0-9]*"),
        self_report_6153 = parse_code_vector(row[["self_report_6153"]], wildcard = "[0-9]*"),
        self_report_6177 = parse_code_vector(row[["self_report_6177"]], wildcard = "[0-9]*"),
        procedures_opcs3 = parse_code_vector(row[["procedures_opcs3"]], wildcard = "[0-9]*"),
        procedures_opcs4 = parse_code_vector(row[["procedures_opcs4"]], wildcard = "[0-9]*"),
        cancer_registry = parse_code_vector(row[["cancer_histology"]], wildcard = "[0-9]*"),
        read_v2 = parse_code_vector(row[["read_v2_primary_care_codes"]], wildcard = "[0-9]*"),
        read_v3 = parse_code_vector(row[["read_v3_primary_care_codes"]], wildcard = "[0-9]*")
      )
    }
    # 2. Apply the function to each row of the raw data, creating a list of lists
    defs_list <- lapply(seq_len(nrow(defs_raw)), function(i) make_def(defs_raw[i, ]))
    
    return(list(defs_raw, defs_list))
  }

# Convenience helper returning a fully parsed configuration object from a
# control sheet path.  The returned list contains:
#   * filename - name for file to create
#   * matched_outcomes – list of matched outcome
#   * diabetes_flag – TRUE if diabetes specific post-processing is required
#   * outcomes_def - list of lists of the definitions defined
parse_control_sheet <- function(path, CVD=F, CANCER=F) {
    settings <- read_settings_sheet(path,CVD = CVD, CANCER = CANCER)
    outcomes <- read_outcome_definitions(path, outcome_filter = settings$outcomes)
    list(
       filename_out = settings$filename,
       matched_outcomes = outcomes$defn_raw,
       diabetes_flag = settings$diabetes_flag,
       outcomes_def = outcomes$defn_list
        )
}
    
    
    
    