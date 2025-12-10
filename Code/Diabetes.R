##########################################
# Use Scott's implementation of diabetes
##########################################

match_diabetes<- function(definitions,
                          pre="Inputs/input_data/prevalent_diabetes.csv",
                          inc="Inputs/input_data/incident_diabetes.csv") {
  library(data.table)
  library(purrr)
  
  
  raw1 <- fread(pre, na.strings = c("", "NA")) 
  
  #this gives me people classified as Diabetes unlikely, Possible gestational diabetes,  
  # Possible type 1 diabetes, Possible type 2 diabetes, Probable type 1 diabetes, Probable type 2 diabetes
  raw2 <- fread(inc, na.strings = c("", "NA")) 
  
  #this gives me people classified as Incident diabetes of uncertain type, No evidence of diabetes,
  # Possible incident type 1 diabetes, Possible incident type 2 diabetes, Prevalent diabetes,
  # Probable incident type 1 diabetes, Probable incident type 2 diabetes
  
  combine <- merge(x=raw1 %>% select(eid, visit_index, adjudicated_diabetes), 
                   y=raw2 %>% select(eid, visit_index, adjudicated_diabetes), 
                   by=c("eid", "visit_index"))
  
  
  
  # consider HES data
  
  # make list of HES codes
  
  outcomes_def= list(
    list(outcome_id="any", icd9="250.X", icd10=c("E10.X", "E11.X", "O24.429")),
    list(outcome_id="gest", icd9="", icd10=c("O24.429")),
    list(outcome_id="T2D", icd9="", icd10=c("E11.X")),
    list(outcome_id="T1D", icd9="", icd10=c("E10.X"))
  )
  
  #hes
  hes_icd9_appears<-load_from_rap(rap_path="/common/Hospital Records/icd9_codes.csv",
                                  work_dir="./Inputs/hes/") %>% #
    filter(appears_in_records==TRUE)%>% pull(code)
  
  hes_icd10_appears<-load_from_rap(rap_path="/common/Hospital Records/icd10_codes.csv",
                                   work_dir="./Inputs/hes/") %>% #
    filter(appears_in_records==TRUE)%>% pull(code)
  
  make_def_diabetes <- function(row) {
    list(
      outcome_id = trimws(as.character(row[["suggested_variable_name"]])),
      icd9 = resolve_code_patterns(code_string = row[["icd_9_codes"]],
                                   catalog_codes = hes_icd9_appears, 
                                   wildcard="[0-9]*"), 
      icd10 = resolve_code_patterns(code_string = row[["icd_10_codes"]],
                                    catalog_codes = hes_icd10_appears, 
                                    wildcard="[0-9]*"), 
      death = resolve_code_patterns(code_string = row[["death_40001_40002"]],
                                    catalog_codes = death_icd10_appears, 
                                    wildcard="[0-9]*"), 
  }
  
  
  defs_list <- lapply(seq_len(nrow(defs_raw)), function(i) make_def(defs_raw[i, ]))
  
  return(list(defs_raw=defs_raw$suggested_variable_name, defs_list=defs_list))
  
  ##
  
  
  
  
  ### load HES, death 
  
  work_dir = "./Inputs"
  
  HES_map = list(
    icd9 = "common/Hospital Records/diagnoses.csv",
    icd10 = "common/Hospital Records/diagnoses.csv",
    death = "common/Deaths/death_causes.csv"  )
  
  rap_paths <- unique(unlist(dataset_map, use.names = FALSE))
  rap_paths <- rap_paths[!is.na(rap_paths) & rap_paths != ""]
  
  loaded <- lapply(rap_paths, load_from_rap, work_dir = work_dir)
  names(loaded) <- basename(rap_paths)
  
  ## create list of needed variables  
  
  load_from_rap
  
  **Add diabetes info from HES
  
  * Definitions for purpose of HES/Death
  
  
  
  
  * SETUP
  **************************************************
    
    * want 4 variables: HES_any, HES_T2, HES_T1, HES_gest
  
  ********************************************
    *HES data
  
  import delimited `HES_diag', clear
capture confirm variable n_eid
if !_rc rename eid n_eid


* match to ICD 10 values

foreach out of newlist gest T1 T2 any{
gen HESmatch_`out' =0
}

* match to ICD 10 values
replace HESmatch_any =1 if regexm(diag_icd10,"E10[0-9]*")
replace HESmatch_any =1 if regexm(diag_icd10,"E11[0-9]*")
replace HESmatch_any =1 if regexm(diag_icd10,"O244")
replace HESmatch_gest =1 if regexm(diag_icd10,"O244")
replace HESmatch_T1 =1 if regexm(diag_icd10,"E10[0-9]*")
replace HESmatch_T2 =1 if regexm(diag_icd10,"E11[0-9]*")

* keep 1 record per person			
keep n_eid HESmatch*
  sort n_eid 
foreach out of newlist any gest T1 T2{
  by n_eid: egen HES_`out' = max (HESmatch_`out')
}
drop HESmatch*
rename HES* HESmatch*
by n_eid: drop if _n>1


*************************************STEP 3: death fields data set
capture confirm variable _merge
if !_rc drop _merge
merge 1:1 n_eid using `inputfile', update
assert inlist(_merge,1,2,3)
drop _merge

tempfile working
save `working', replace

import delimited `DEATH_add_date', clear
capture confirm variable eid
if !_rc rename n_eid eid

* reshape to wide: 1 record per eid
keep eid ins date
reshape wide date, i(eid) j(ins)

* check second records of dates don't introduce multiple options for death date
assert date_of_death0 ==date_of_death1 if date_of_death1!=""

* turn into date format
gen ts_40000_2_0 = date(date_of_death0, "DMY")
format ts* %td
drop date*
*   
* save in temp file
tempfile deathtemp
save `deathtemp', replace

** load second dataset
import delimited `DEATH_add_cause', clear

* reshape to wide
keep eid level cause
rename cause cause
sort eid level cause
bysort eid level: gen count=_n
reshape wide cause, i(eid level) j(count)
rename cause* cause*_
sort eid level
reshape wide cause*, i(eid) j(level)
rename cause*_1 s_40001_2_*
rename cause*_2 s_40002_2_*


* merge with data dataset
merge 1:1 eid using `deathtemp', update
drop _merge

keep eid HES* s_40001* s_40002*
  
  ****** death report 40001

gen Deathmatch_any =0
gen Deathmatch_gest =0
gen Deathmatch_T1=0
gen Deathmatch_T2=0

*match to ICD 10 values
foreach i of varlist s_40001* s_40002*{
  replace Deathmatch_any =1 if regexm(`i',"E10[0-9]*")
replace Deathmatch_any =1 if regexm(`i',"E11[0-9]*")
  replace Deathmatch_any =1 if regexm(`i',"O244")
replace Deathmatch_gest =1 if regexm(`i',"O244")
  replace Deathmatch_T1 =1 if regexm(`i',"E10[0-9]*")
replace Deathmatch_T2 =1 if regexm(`i',"E11[0-9]*")
}


drop s*
  
  rename eid n_eid

capture confirm variable _merge
if !_rc drop _merge

merge 1:1 n_eid using `working', update
drop _merge


* combine death and HES
gen HESDeath_dia_any = Deathmatch_any
replace HESDeath_dia_any = 1 if HESmatch_any==1

gen HESDeath_dia_gest = Deathmatch_gest
replace HESDeath_dia_gest = 1 if HESmatch_gest==1

gen HESDeath_dia_T1 = Deathmatch_T1
replace HESDeath_dia_T1 = 1 if HESmatch_T1==1

gen HESDeath_dia_T2 = Deathmatch_T2
replace HESDeath_dia_T2 = 1 if HESmatch_T2==1

* create single variable

gen summary_HESDeath = 1 if HESDeath_dia_any==0
replace summary_HESDeath = 0 if HESDeath_dia_any==1
replace summary_HESDeath = 2 if HESDeath_dia_gest == 1 & HESDeath_dia_T1 ==0 & HESDeath_dia_T2 ==0 
replace summary_HESDeath = 6 if HESDeath_dia_T1 ==1 & HESDeath_dia_T2 ==0 
replace summary_HESDeath = 4 if HESDeath_dia_T1 ==0 & HESDeath_dia_T2 ==1 
lab val summary_H summary_def

noi di "HES and DEath info"
noi tab summary

keep n_eid summary_HESDeath
  
  # I wnat people as T2D_probable : yes/no/NA
  # T2D_possible : yes/no/NA
  # T1D_probable:
  # T1D_possible: 
  
  