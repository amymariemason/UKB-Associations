## function to match self-report

## input outcome definitions
## location of self_report file



match_self_report<-function(definitions=test2$outcomes_def, 
                            self_report_file="~/Inputs/ukb_input_data/data.csv", 
                            suffix="_SR"){
  
  # Load in data
  library(data.table)
  library(purrr)
  
  raw <- fread(self_report_file, na.strings=c("", "NA"))
  
  # list of self report fields
  
  # Manually curate field information we want to extract
  info <- rbind(use.names=FALSE,
                data.table(field.id=2966, name="Age high blood pressure diagnosed (touchscreen)"),
                data.table(field.id=2976, name="Age diabetes diagnosed (touchscreen)"),
                data.table(field.id=2986, name="Started insulin within one year diagnosis of diabetes (touchscreen)"),
                data.table(field.id=3627, name="Age angina diagnosed (touchscreen)"),
                data.table(field.id=3761, name="Age hay fever, rhinitis or eczema diagnosed (touchscreen)"),
                data.table(field.id=3786, name="Age asthma diagnosed (touchscreen)"),
                data.table(field.id=3894, name="Age heart attack diagnosed (touchscreen)"),
                data.table(field.id=3992, name="Age emphysema/chronic bronchitis diagnosed (touchscreen)"),
                data.table(field.id=4012, name="Age deep-vein thrombosis (DVT, blood clot in leg) diagnosed (touchscreen)"),
                data.table(field.id=4022, name="Age pulmonary embolism (blood clot in lung) diagnosed (touchscreen)"),
                data.table(field.id=4041, name="Gestational diabetes only (touchscreen)"),
                data.table(field.id=4056, name="Age stroke diagnosed (touchscreen)"),
                data.table(field.id=6150, name="Vascular/heart problems diagnosed by doctor (touchscreen)"),
                data.table(field.id=6152, name="Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor (touchscreen)"),
                data.table(field.id=6177, name="Touchscreen medication question for males"),
                data.table(field.id=6153, name="Touchscreen medication question for females"),
                data.table(field.id=20001, name="Cancer code (verbal interview)"),
                data.table(field.id=20002, name="Non-cancer illness code (verbal interview)"),
                data.table(field.id=20003, name="Verbal interview with nurse: treatment/medication code"),
                data.table(field.id=20004, name="Operation code (verbal interview)"),
                data.table(field.id=20006, name="Interpolated year when cancer first diagnosed (verbal interview)"),
                data.table(field.id=20007, name="Interpolated age when cancer first diagnosed (verbal interview)"),
                data.table(field.id=20008, name="Interpolated year when non-cancer illness first diagnosed (verbal interview)"),
                data.table(field.id=20009, name="Interpolated age when non-cancer illness first diagnosed (verbal interview)"),
                data.table(field.id=20010, name="Interpolated year when operation took place (verbal interview)"),
                data.table(field.id=20011, name="Interpolated age when operation took place (verbal interview)")
  )
  
  
  # break down to list per outcome
  
  # Convert to long format per field - we split this out into a list, one per 
  # field, as fields are a mix of codes (integer), age of event (numeric), and
  # date of event.
  self_report <- lapply(info$field.id, function(fid) {
    dt <- raw[,.SD,.SDcols=c("eid", names(raw)[names(raw) %like% sprintf("^p%s_", fid)])]
    suppressWarnings(dt <- melt(dt, id.vars="eid", na.rm=TRUE)) # don't warn on coercion
    field_visit_repeat <- dt[,tstrsplit(variable, "_")]
    if (ncol(field_visit_repeat) == 2) field_visit_repeat[, V3 := "a0"]
    setnames(field_visit_repeat, c("field_id", "visit_index", "repeat_index"))
    dt <- cbind(dt, field_visit_repeat)
    if (is.character(dt$value)) {
      # sometimes multiple choice answers are collapsed into a single entry - make multiple rows
      dt <- dt[, .(value=strsplit(value, "\\|")[[1]]), by=.(eid, visit_index, repeat_index)]
    }
    dt[, visit_index := as.integer(gsub("i", "", visit_index))]
    dt[, repeat_index := as.integer(gsub("a", "", repeat_index))]
    dt <- dt[, .(eid, visit_index, repeat_index, value=as.numeric(value))]
    return(dt)
  })
  names(self_report) <- as.character(info$field.id)
  
  # combine the self report fields only
  sr <- rbindlist(self_report, use.names=TRUE, idcol="field_id")
  sr$field_id<-as.numeric(sr$field_id)
  
  
  # flatten the self-report definitions
  deflist <- rbindlist(
    lapply(definitions, function(def) {
      
      # find elements with names like self_report_XXXX
      sr_fields <- names(def)[grepl("^self_report_", names(def))]
      
      if (length(sr_fields) == 0) return(NULL)
      
      rbindlist(lapply(sr_fields, function(nm) {
        field_id <- as.integer(sub("self_report_", "", nm))
        
        # patterns may be character(0) or NULL
        patterns <- def[[nm]]
        if (length(patterns) == 0) return(NULL)
        
        data.table(
          outcome_id = def$outcome_id,
          field_id   = field_id,
          pattern    = patterns
        )
      }))
    }),
    use.names=TRUE
  )
  
  if (nrow(deflist)==0) {
    stop("No self-report definitions found in the provided definition set.")
  }
  
  
  # filter to only the requested SR ids
  
  unique_sr<- unique(deflist$field_id)
  sr<- sr %>% filter(field_id %in% unique_sr)
  
  
  # join the sets
  
  setkey(sr, field_id)
  setkey(deflist, field_id)
  
  merged <- deflist[sr, nomatch=0L, allow.cartesian=TRUE]
  
  ## check for matches
  
  merged[, match := mapply(grepl, pattern, value)]
  
  
  ## collapse to binary indicator for each eid and outcome
  ## NB: this only contains the matches
  indicator <- merged[
    match == TRUE,
    .(has_self_report = TRUE),
    by = .(eid, outcome_id)
  ]
  
  rm(merged)
  
  ## fill in all the outcomes
  
  wide <- dcast(
    indicator,
    eid ~ outcome_id,
    value.var = "has_self_report",
    fill = FALSE
  )
  
  ### fill in all eids
  
  all_eids <- unique(raw$eid)
  
  
  final <- merge(
    data.table(eid = all_eids),
    wide,
    by = "eid",
    all.x = TRUE
  )
  
assert(length(final$eid)==length(raw$eid))

# fill back with false where no report found

  for (col in names(final)[-1]) {
    set(final, which(is.na(final[[col]])), col, FALSE)
  }

  ## add suffix to outcome
  if (suffix != "") {
    newnames <- paste0(names(final)[-1], suffix)
    setnames(final, old = names(final)[-1], new = newnames)
  }
  
  # ---- 9. On-screen summary ----
  message("\nSelf-report match summary:")
  counts <- final[, lapply(.SD, sum), .SDcols = names(final)[-1]]
  print(counts)

invisible(final)
}

