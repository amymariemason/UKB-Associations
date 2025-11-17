#### This is based on some helpful code from scott richie
#
#
# Function to load a dataset from the RAP
load_from_rap <- function(rap_path, work_dir) {
  if (!dir.exists(sprintf("%s/input_data/", work_dir))) {
    system(sprintf("mkdir -p %s/input_data", work_dir), ignore.stdout=FALSE)
  }
  fname <- basename(rap_path)
  if (!file.exists(sprintf("%s/input_data/%s", work_dir, fname))) {
    system(sprintf("dx download '%s' -o '%s/input_data/%s'", rap_path, work_dir, fname), ignore.stdout=FALSE)
  }
  fread(sprintf("%s/input_data/%s", work_dir, fname), na.strings=c("", "NA"))
}

# Load ICD codes associated with cancer register
cancer_icd10 <- load_from_rap("common/Cancer Register/icd10_codes.csv","./Inputs")
cancer_icd9 <- load_from_rap("common/Cancer Register/icd9_codes.csv")