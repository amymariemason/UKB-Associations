### working flowchart

# source code
source("~/Code/Read_settings.R")
source("~/Code/gather_UKB_files.R")

# input definitions location
path<- "Inputs/bespoke_outcome_v4.xlsx"  

#make definitions long file
test<- parse_control_sheet(path, CVD=T)   

# create summary of what data is needed
test2<- curate_settings(test)

outcome_list_by_dataset <- test2$data_loop_lists
datasets_needed<-test2$data_requirements

load_required_datasets(test2)

# create self-report variables

source("~/Code/self_report.R")
outcomes_SR<-match_self_report(definitions=test2$outcomes_def, 
                            self_report_file="~/Inputs/ukb_input_data/data.csv", 
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



