################################################################################
#               Definitions Pipeline for MR
################################################################################


# list outcomes wanted
# this should be "CVD" or "Cancer" or NULL or a vector of characters

CVD <- TRUE
CANCER <- FALSE
SHEET <-FALSE  # take the ticked list from an edited sheet
INPUT <-  NULL #  # list as vector e.g. c("ca_all", "ca_hep", "aa")
DIABETES<-FALSE


custom_filename= "CVD"


# input definitions location
path<- "Inputs/bespoke_outcome_v4.xlsx"  

###### Source code

# source code
source("~/Code/Read_settings.R")
source("~/Code/gather_UKB_files.R")

###### Collect definitions


    #make definitions long file
    def_list<- parse_control_sheet(path, 
                               CANCER=CANCER, 
                               CVD=CVD, 
                               SHEET=SHEET,
                               INPUT=INPUT,
                               DIABETES= DIABETES,
                               custom_filename= custom_filename)   
    
    # create summary of what data is needed
    settings<- curate_settings(def_list)
    
    outcome_list_by_dataset <- settings$data_loop_lists
    datasets_needed<-settings$data_requirements
    
    load_required_datasets(settings)
    
    outcomes_def <- settings$outcomes_def

###### Collect data from sources
    
    # create initial list of eids
    
    eid_list<-load_from_rap("/users/Amy/MR_base_european_unrelated.csv",work_dir = "./Inputs")
    
    # create self-report variables
    if (datasets_needed$self_report) {
    source("~/Code/self_report.R")
    outcomes_SR<-match_self_report(definitions=outcomes_def, 
                                   self_report_file="~/Inputs/input_data/data.csv", 
                                   suffix="_SR")
    }
    
    ## create hes variables
    if (datasets_needed$icd9 |datasets_needed$icd10) {
    source("~/Code/hes.R")
    outcomes_hes<- match_hes1(definitions=outcomes_def, 
                              hospital_file="Inputs/input_data/diagnoses.csv",
                              suffix="_hes")
    
    }
    
    if (datasets_needed$procedures){
    source("~/Code/hes.R")
    outcomes_opcs<- match_hes2(definitions=outcomes_def, 
                               hospital_file="Inputs/input_data/operations.csv",
                               suffix="_opcs")
    }
    ## create death variables
    
    if (datasets_needed$death){
    source("~/Code/death.R")
    outcomes_death <- match_death(definitions=outcomes_def, 
                                  death_file="Inputs/input_data/death_causes.csv",
                                  suffix="_death")
    }
    ## create cancer variables
    if (datasets_needed$cancer_registry){
      source("~/Code/cancer.R")
      
      if (datasets_needed$icd10){
      outcomes_cancer <- match_cancer_icd(definitions=outcomes_def, 
                                    cancer_file="Inputs/input_data/cancer_register.csv",
                                    suffix="_cancer_icd") 
      }
      
      outcomes_histology<- match_cancer_histology(definitions=outcomes_def, 
                                    cancer_file=  "Inputs/input_data/cancer_register.csv",
                                    suffix="_cancer_hist") 
      
    }
    
    ## create primary care 
    if (datasets_needed$primary_care_read){
      stop("ERROR: primary care not yet implemented")
      outcomes_primary<-NULL
    }
    
    
###### Merge all created datasets
    
    # List of potential dataframe names created by the sourcing steps
    possible_dfs <- c("outcomes_SR", "outcomes_hes", "outcomes_opcs", "outcomes_death", "outcomes_cancer", "outcomes_histology", "outcomes_primary")
    
    # Filter to keep only the dataframes that actually exist in the environment
    existing_dfs <- possible_dfs[sapply(possible_dfs, exists)]
    
    # Merge them all by 'eid'
    # We use Reduce() to merge a list of dataframes sequentially
    if(length(existing_dfs) > 0) {
      
      # Put actual dataframes into a list
      list_of_dfs <- mapply(get, existing_dfs, SIMPLIFY = FALSE)
      
      # Merge all existing dataframes by 'eid' (all=TRUE ensures we keep all people even if they don't have data in one source)
      merged_outcomes <- Reduce(function(x, y) merge(x, y, by = "eid", all = TRUE), list_of_dfs)
      
      
    } else {
      stop("No outcome dataframes were created.")
    }

###### Create Summary Variables
    
    # Get unique disease names from your definitions
    disease_names <- def_list$matched_outcomes
    
    # Loop through each disease to create the summary columns
    for (disease in disease_names) {
      
      # 1. Identify all columns related to this disease
      cols_all <- grep(paste0("^", disease), names(merged_outcomes), value = TRUE)
      
      # 2. Identify columns EXCLUDING self-report
      cols_no_sr <- cols_all[!grepl("_SR$", cols_all)]
      
      # 3. Create "Any Source" variable
      if(length(cols_all) > 0){
        merged_outcomes[[paste0(disease, "_all_reports")]] <- rowSums(merged_outcomes[, ..cols_all, drop=FALSE], na.rm = TRUE) > 0
      } else {
        merged_outcomes[[paste0(disease, "_all_reports")]] <- FALSE
      }
      
      # 4. Create "Any Source Excluding Self-Report" variable
      if(length(cols_no_sr) > 0){
        merged_outcomes[[paste0(disease, "_nSR")]] <- rowSums(merged_outcomes[, ..cols_no_sr, drop=FALSE], na.rm = TRUE) > 0
      } else {
        # If there are no clinical columns (e.g. only SR existed), this is FALSE
        merged_outcomes[[paste0(disease, "_nSR")]] <- FALSE
      }
    }
    
    # keep only the all and nSR variables
    merged_outcomes<- merged_outcomes %>%
      select(eid, sort(ends_with("_all_reports")), sort(ends_with("_nSR"))) %>%
      relocate(eid) %>%
      rename_with(~ str_remove(., "_all_reports$"), ends_with("_all_reports"))

    
# save in outputs
    
#     

 write_csv(merged_outcomes, file=paste0("./Outputs/",settings$filename_out ,".csv"))
 
# transfer back to RAP
 
 system(sprintf("dx upload '%s' -o './users/Amy/Outcome files/%s'", 
                paste0("./Outputs/",settings$filename_out ,".csv"), 
                paste0(settings$filename_out ,".csv")), ignore.stdout = FALSE)
 