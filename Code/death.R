###############
### match death

### as with the self report, but this time use Scott's preprocessed data

match_death <- function(definitions,
                     death_file="Inputs/input_data/death_causes.csv",
                     suffix="_death") {
  library(data.table)
  library(purrr)
  
  
  raw <- fread(death_file, na.strings = c("", "NA")) 
  
 death<- raw %>%
    select(eid, cause_icd10)
  all_eids <-unique(raw$eid)
  
  # make definition list
  deflist <- rbindlist(
    lapply(definitions, function(def) {
      icd_fields <- names(def)[grepl("^death", names(def))]
      if (!length(icd_fields)) return(NULL)
      
      rbindlist(lapply(icd_fields, function(nm) {
        field_icd_type <- nm
        pats <- def[[nm]]
        if (is.null(pats) || length(pats) == 0) return(NULL)
        
        data.table(
          outcome_id = def$outcome_id,
          cause_icd10   = as.character(pats)
        )
      }), use.names = TRUE)
    }),
    use.names = TRUE, fill = TRUE
  )
  
  
  if (nrow(deflist) == 0) {
    stop("No hes patterns found in definitions.")
  }
  
  results <- list()
  
    message("Matching ICD10 to death certificates ")
    
    # Exact match via join (FAST!)
    setkey(death, cause_icd10)
    setkey(deflist, cause_icd10)
    
    merged <- deflist[death, nomatch = 0L, allow.cartesian = TRUE]
    
    if (nrow(merged) == 0) {
      message("No deaths found")
      return(data.table(eid = all_eids))
    }
    
    # Collapse to binary indicator
    results <- merged[, .(has = TRUE), by = .(eid, outcome_id)]
    

  # ---- Wide reshape ----
  wide <- dcast(results, eid ~ outcome_id, value.var = "has", fill = FALSE)
  
  # ---- Include all participants ----
  final <- merge(data.table(eid = all_eids), wide, by="eid", all.x=TRUE)
  
  for (col in names(final)[-1]) {
    set(final, which(is.na(final[[col]])), col, FALSE)
  }
  
  # ---- Add _SR suffix ----
  if (suffix != "") {
    setnames(final, old = names(final)[-1],
             new = paste0(names(final)[-1], suffix))
  }
  # ---- Diagnostics ----
  message("\nICD match counts:")
  print(final[, lapply(.SD, sum), .SDcols = names(final)[-1]])
  
  return(final)  
  
}



