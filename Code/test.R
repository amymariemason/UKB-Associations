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
