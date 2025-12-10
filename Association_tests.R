# === R script to run PLINK2 logistic per-phenotype and import results ===
library(data.table)   # fread/fwrite speed
library(dplyr)

# 1) paths - edit these
bfile_prefix <- "alcohol_and_smoking"   # points to mydata.bed/mydata.bim/mydata.fam
pheno_csv    <- "phenotypes/phenos.csv"  # your CSV: must contain columns for FID/IID or sample ID
plink_exec   <- "~/_ukbrapr_tolls/plink2"  # or full path, e.g. "/usr/local/bin/plink2"
out_dir      <- "plink_results"
dir.create(out_dir, showWarnings = FALSE)

# 2) read fam and your phenotype CSV
fam <- fread(paste0(bfile_prefix, ".fam"), header = FALSE)
# .fam columns: FID IID PID MID SEX PHEN (but headerless)
colnames(fam)[1:2] <- c("FID", "IID")

phenos <- fread(pheno_csv)   # adjust below mapping to your CSV column names

# ---- Example: if your CSV has columns "eid" matching IID, or columns FID/IID ----
# If your CSV has only one ID column (IID), make FID=IID:
if(!("FID" %in% names(phenos)) && ("IID" %in% names(phenos))){
  phenos <- phenos %>% mutate(FID = IID)
}
if(!("FID" %in% names(phenos)) && ("eid" %in% names(phenos))){
  phenos <- phenos %>% rename(IID = eid) %>% mutate(FID = IID)
}

# 3) merge to ensure consistent samples (keeps order by FID/IID to match fam if needed)
merged <- fam %>% select(FID, IID) %>% left_join(phenos, by = c("FID", "IID"))

# 4) write a PLINK-style phenotype file (header optional; PLINK2 accepts header)
#    We'll write FID, IID, then multiple phenotype columns named P1, P2 ... OR keep original names
# Suppose your case-control columns are named pheno1, pheno2, etc.; adjust accordingly.
case_cols <- grep("pheno", names(merged), value = TRUE) # adjust pattern
if(length(case_cols) == 0) stop("No phenotype columns found: adjust `case_cols` selection")

pheno_out <- merged %>% select(FID, IID, all_of(case_cols))
fwrite(pheno_out, file = file.path(out_dir, "plink_pheno.txt"), sep = "\t", na = "NA")

# 5) write covariates file (example: age, sex, PC1..PC10)
covar_cols <- c("age", "sex", paste0("PC", 1:10))  # edit based on your data
covar_present <- intersect(covar_cols, names(merged))
if(length(covar_present) > 0){
  covar_out <- merged %>% select(FID, IID, all_of(covar_present))
  fwrite(covar_out, file = file.path(out_dir, "covar.txt"), sep = "\t", na = "NA")
} else {
  message("No covariates found; proceeding without covariates.")
}

# 6) function to call plink2 for one phenotype name
run_plink_for_pheno <- function(pheno_name){
  out_prefix <- file.path(out_dir, paste0("plink_", pheno_name))
  cmd <- sprintf("%s --bfile %s --pheno %s --pheno-name %s",
                 plink_exec, bfile_prefix,
                 file.path(out_dir, "plink_pheno.txt"),
                 pheno_name)
  # add covar if present
  if(length(covar_present) > 0){
    cmd <- paste(cmd, sprintf("--covar %s --covar-name %s",
                              file.path(out_dir, "covar.txt"),
                              paste(covar_present, collapse=",")))
  }
  # association: glm with firth fallback, and request relevant cols (odds ratio etc)
  cmd <- paste(cmd, "--glm firth-fallback cols=+a1freq,odds_ratio --allow-no-sex --out", out_prefix)
  message("Running: ", cmd)
  status <- system(cmd)
  if(status != 0) warning("plink returned non-zero status for ", pheno_name)
  # Read sumstats: PLINK2 default .glm.* or .sscore / .ssv - choose the file produced
  # PLINK2 produces <out_prefix>.glm.logistic (or .glm) - inspect actual output name for your version.
  glmfile <- paste0(out_prefix, ".glm.logistic") 
  if(file.exists(glmfile)){
    res <- fread(glmfile, na.strings = c("NA","nan"))
    res$PHENO <- pheno_name
    return(res)
  } else {
    # Try GWAS-SSF output if configured
    ssf <- paste0(out_prefix, ".ssf.tsv")
    if(file.exists(ssf)) {
      res2 <- fread(ssf)
      res2$PHENO <- pheno_name
      return(res2)
    } else {
      warning("No GLM output found for ", pheno_name)
      return(NULL)
    }
  }
}

# 7) loop through phenotypes and collect results
all_results <- list()
for(ph in case_cols){
  res <- run_plink_for_pheno(ph)
  if(!is.null(res)) all_results[[ph]] <- res
}

# 8) combine and save
if(length(all_results) > 0){
  combined <- rbindlist(all_results, fill = TRUE)
  fwrite(combined, file = file.path(out_dir, "combined_results.tsv"), sep = "\t")
  message("Results combined: ", file.path(out_dir, "combined_results.tsv"))
} else {
  message("No results to combine.")
}
