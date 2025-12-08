################################################################################
#               Definitions Pipeline for MR
################################################################################


# list outcomes wanted
# this should be "CVD" or "Cancer" or NULL or a vector of characters

CVD <- TRUE
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
    test<- parse_control_sheet(path, 
                               CANCER=CANCER, 
                               CVD=CVD, 
                               SHEET=SHEET,
                               INPUT=INPUT,
                               custom_filename= custom_filename)   
    
    # create summary of what data is needed
    test2<- curate_settings(test)
    
    outcome_list_by_dataset <- test2$data_loop_lists
    datasets_needed<-test2$data_requirements
    
    load_required_datasets(test2)

###### Collect data from sources
    
    # create initial list of eids
    
    eid_list<-load_from_rap("/users/Amy/MR_base_european_unrelated.csv",work_dir = "./Inputs")
    
    # create self-report variables
    
    source("~/Code/self_report.R")
    outcomes_SR<-match_self_report(definitions=test2$outcomes_def, 
                                   self_report_file="~/Inputs/input_data/data.csv", 
                                   suffix="_SR")
    
    ## create hes variables
    source("~/Code/hes.R")
    outcomes_hes<- match_hes1(definitions=test2$outcomes_def, 
                              hospital_file="Inputs/input_data/diagnoses.csv",
                              suffix="_hes")
    
    outcomes_opcs<- match_hes2(definitions=test2$outcomes_def, 
                               hospital_file="Inputs/input_data/operations.csv",
                               suffix="_opcs")
    
    ## create death variables
    
    source("~/Code/death.R")
    outcomes_death <- match_death(definitions=test2$outcomes_def, 
                                  death_file="Inputs/input_data/death_causes.csv",
                                  suffix="_death")
    
    ## create cancer variables
    
    
    
    ## create primary care 


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