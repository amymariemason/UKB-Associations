###############
### match cancer

### as with the self report, but this time use Scott's preprocessed data

match_cancer_icd<- function(definitions,
                        cancer_file="Inputs/input_data/cancer_register.csv",
                        suffix="_cancer") {
  library(data.table)
  library(purrr)
  
  
  raw <- fread(cancer_file, na.strings = c("", "NA")) 

  cancer<- raw %>%
    select(eid, icd_version, diag_icd)
  all_eids <-unique(raw$eid)
  
  # make definition list
  deflist <- rbindlist(
    lapply(definitions, function(def) {
      icd_fields <- names(def)[grepl("^icd", names(def))]
      if (!length(icd_fields)) return(NULL)
      
      rbindlist(lapply(icd_fields, function(nm) {
        field_icd_type <- sub("icd","",nm)
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
    stop("No icd patterns found in definitions.")
  }
  
  # Ensure types
  cancer[, icd_version := as.integer(icd_version)]
  deflist[, icd_version := as.integer(icd_version)]
  
  
  # ---- Split by ICD version (ICD-9 vs ICD-10) ----
  versions <- sort(unique(deflist$icd_version))
  
  results <- list()
  
  for (v in versions) {
    message("Matching Cancer Register ICD version ", v, "...")
    
    def_v <- deflist[icd_version == v]
    cancer_v <- cancer[icd_version == v]
    
    if (nrow(cancer_v) == 0 || nrow(def_v) == 0) next
    
    # Exact match via join (FAST!)
    setkey(cancer_v, diag_icd)
    setkey(def_v, pattern)
    
    merged <- def_v[cancer_v, nomatch = 0L, allow.cartesian = TRUE]
    
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
  wide <- dcast(indicator_all, eid ~ outcome_id, 
                value.var = "has", 
                fill = FALSE,
                fun.aggregate = any)
  
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
  message("\nCancer register ICD code match counts:")
  print(final[, lapply(.SD, sum), .SDcols = names(final)[-1]])
  
  return(final)  
  
}

match_cancer_histology<- function(definitions,
                            cancer_file="Inputs/input_data/cancer_register.csv",
                            suffix="_cancer_hist") {
  library(data.table)
  library(purrr)
  
  
  raw <- fread(cancer_file, na.strings = c("", "NA")) 
  
  cancer<- raw %>%
    mutate(hist=sub("/.*", "",sub("M-", "", icdO3_code)),
           type= "hist_code") %>%
    select(eid, hist, type)%>%
    filter(!(is.na(hist)|hist==""))
    
  all_eids <-unique(raw$eid)
  
  # make definition list
  deflist <- rbindlist(
    lapply(definitions, function(def) {
      hist_fields <- names(def)[grepl("^cancer", names(def))]
      if (!length(hist_fields)) return(NULL)
      
      rbindlist(lapply(hist_fields, function(nm) {
        pats <- def[[nm]]
        if (is.null(pats) || length(pats) == 0) return(NULL)
        
        data.table(
          outcome_id = def$outcome_id,
          type="hist_code",
          pattern    = as.character(pats)
        )
      }), use.names = TRUE)
    }),
    use.names = TRUE, fill = TRUE
  )
  
  
  if (nrow(deflist) == 0) {
    stop("No histology patterns found in definitions.")
  }
  
  
  # join the sets
  
  setkey(cancer, type)
  setkey(deflist, type)
  
  merged <- deflist[cancer, nomatch=0L, allow.cartesian=TRUE]
  
  ## check for matches
  
  merged[, match := mapply(grepl, pattern, hist)]
  
  
  ## collapse to binary indicator for each eid and outcome
  ## NB: this only contains the matches
  indicator <- merged[
    match == TRUE,
    .(has = TRUE),
    by = .(eid, outcome_id)
  ]
 
  
  # ---- Combine all ICD versions ----
  if (length(results) == 0) {
    warning("No histology matches found.")
    return(data.table(eid = all_eids))
  }
  
  # ---- Wide reshape ----
  wide <- dcast(indicator, eid ~ outcome_id, 
                value.var = "has", 
                fill = FALSE,
                fun.aggregate = any)
  
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
  message("\nCancer register histology match counts:")
  print(final[, lapply(.SD, sum), .SDcols = names(final)[-1]])
  
  return(final)  
  
}

