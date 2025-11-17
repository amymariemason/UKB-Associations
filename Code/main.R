
# TO DO - use curated icd10 files to keep exact match options only in Read_settings
# find all the relevant files and add them to path

# will need an extra version of this that always loads the icd curation lists

# loop over files and create matches


### load needed functions
source("~/UKB-Associations/Code/gather_UKB_files.R")
source("~/UKB-Associations/Code/Read_settings.R")

## curate list of required functions and data
path<- "Inputs/bespoke_outcome_v3.xls"
required<-curate_settings(parse_control_sheet(path, CVD=T, CANCER = F))

## gather needed files
rap_paths <- list(
  icd9 = "path/to/icd9.csv",
  icd10 = "path/to/icd10.csv",
  death = "path/to/death.csv",
  self_report = "path/to/self_report.csv",
  procedures = "path/to/opcs.csv",
  cancer_registry = "CEU_overarching/common/Cancer Register/cancer_register.csv",
  primary_care_read = "path/to/gpread.csv"
)