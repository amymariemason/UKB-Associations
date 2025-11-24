# based on Scott's curate field file

library(data.table)
library(foreach)

# Check if we have the necessary python libraries installed for dx extract_dataset
if (system("python3 -c 'import pandas'", ignore.stderr=TRUE)) {
  system("pip install pandas", wait=TRUE)
}

# Extract dataset meta-data to local cloud instance, if it hasn't already been downloaded
if (length(list.files(pattern="data_dictionary")) == 0) {
  project_files <- system("dx ls", intern=TRUE)
  dataset_file <- project_files[project_files %like% "^app" & project_files %like% ".dataset$"]
  system(sprintf("dx extract_dataset %s -ddd", dataset_file), wait=TRUE) # Can ignore warning about pandas version
}

# Auto-detect data-dictionary file then load
data_dict_file <- list.files(pattern="data_dictionary")
data_dict <- fread(data_dict_file)

# Manually curate field information we want to extract
info <- rbind(use.names=FALSE,
              data.table(field.id=31, name="Sex"),
              data.table(field.id=34, name="Year of birth"),
              data.table(field.id=52, name="Month of birth"),
              data.table(field.id=53, name="Date of attending assessment centre"),
              data.table(field.id=54, name="Assessment centre"),
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
              data.table(field.id=20011, name="Interpolated age when operation took place (verbal interview)"),
              data.table(field.id=21000, name="Ethnicity"),
              data.table(field.id=21003, name="Age at assessment"),
              data.table(field.id=22189, name="Townsend deprivation index at recruitment"),
              data.table(field.id=40005, name="Cancer date (cancer registry)"),
              data.table(field.id=40006, name="Cancer type icd10 (cancer registry)"),
              data.table(field.id=40008, name="Age at cancer diagnosis (interpolated)"),
              data.table(field.id=40011, name="Cancer histology (cancer registry)"),
              data.table(field.id=40013, name="Cancer type icd9 (cancer registry)")
              )

# Save table of information
system("mkdir -p medical_history")
fwrite(info, file="medical_history/Amy_field_information.csv")

# Extract list of fields for Table Exporter
fields <- foreach(this_field_id = info$field.id, .combine=rbind) %do% {
  data_dict[name %like% sprintf("p%s_", this_field_id)]
}

# Add in participant ID
fields <- rbind(data_dict[entity == "participant" & name == "eid"], fields)

# Write out field IDs for each distinct entity (only participant in this case,
# but generalizes later if we want to extract different entity types)
for (entity_type in unique(fields$entity)) {
  fwrite(fields[entity == entity_type, .(name)], quote=FALSE, col.names=FALSE,
         file=sprintf("medical_history/Amy_medical_history_fields_%s_entity.txt", entity_type))
}

# Upload field list to persistent storage
system("dx upload medical_history/Amy* --destination 'users/Amy/Define Outcomes Script/'", wait=TRUE)
