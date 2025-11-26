###############
### match death

### as with the self report, but this time use Scott's preprocessed data

make_hes <- function(definitions,
                     death_file="Inputs/input_data/death_causes.csv",
                     suffix="_hes") {
  library(data.table)
  library(purrr)
  
  
  raw <- fread(death_file, na.strings = c("", "NA")) 
  
 death<- raw %>%
    select(eid, icd_version, diag_icd)
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
          icd_version   = field_icd_type,
          pattern    = as.character(pats)
        )
      }), use.names = TRUE)
    }),
    use.names = TRUE, fill = TRUE
  )
  
  
  if (nrow(deflist) == 0) {
    stop("No hes patterns found in definitions.")
  }
  
  # Ensure types
  death[, icd_version := as.integer(icd_version)]
  deflist[, icd_version := as.integer(icd_version)]
  
  
  # ---- Split by ICD version (ICD-9 vs ICD-10) ----
  versions <- sort(unique(deflist$icd_version))
  
  results <- list()
  
  for (v in versions) {
    message("Matching ICD version ", v, "...")
    
    def_v <- deflist[icd_version == v]
    hosp_v <- hosp[icd_version == v]
    
    if (nrow(hosp_v) == 0 || nrow(def_v) == 0) next
    
    # Exact match via join (FAST!)
    setkey(hosp_v, diag_icd)
    setkey(def_v, pattern)
    
    merged <- def_v[hosp_v, nomatch = 0L, allow.cartesian = TRUE]
    
    if (nrow(merged) == 0) next
    
    # Collapse to binary indicator
    indicator <- merged[, .(has = TRUE), by = .(eid, outcome_id)]
    
    results[[as.character(v)]] <- indicator
  }
  
  # ---- Combine all ICD versions ----
  if (length(results) == 0) {
    warning("No ICD matches found.")
    return(data.table(eid = all_eids))
  }
  
  indicator_all <- rbindlist(results, use.names = TRUE)
  
  # ---- Wide reshape ----
  wide <- dcast(indicator, eid ~ outcome_id, value.var = "has", fill = FALSE)
  
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



