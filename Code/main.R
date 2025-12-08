################################################################################
#               Definitions Pipeline for MR
################################################################################


# list outcomes wanted
# this should be "CVD" or "Cancer" or NULL or a vector of characters

CVD <- FALSE
CANCER <- FALSE
SHEET <-FALSE  # take the ticked list from an edited sheet
INPUT <- c("ca_all", "ca_hep", "aa") # list as vector
DIABETES<-TRUE


custom_filename= NULL


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
      
      if (datasets_needed$icd10){
      outcomes_cancer <- match_cancer_icd(definitions=outcomes_def, 
                                    death_file="Inputs/input_data/death_causes.csv",
                                    suffix="_cancer_icd") 
      }
      
      outcomes_histology<- match_cancer_histology(definitions=outcomes_def, 
                                          death_file="Inputs/input_data/death_causes.csv",
                                          suffix="_cancer_hist") 
      
    }
    
    ## create primary care 
    if (datasets_needed$cancer_registry){
      stop("ERROR: primary care not yet implemented")
    }
    
    

