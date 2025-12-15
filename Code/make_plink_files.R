### reshape to plink2 output

## get fam file

fam <- read.table("~/alcohol_and_smoking.fam", quote="\"", comment.char="")

colnames(fam) <- c(
  "FID", "IID", "PID", "MID", "SEX", "PHENO"
)



#### covariates

covar<- load_from_rap("/users/Amy/MR_base_european_unrelated.csv",work_dir = "./Inputs")

covar_plink <- covar %>%   
  transmute(
    FID = as.numeric(eid),
    IID = as.numeric(eid),
    PC1, PC2, PC3, PC4, PC5,
    PC6, PC7, PC8, PC9, PC10,
    sex = case_when(
      sex=="Male" ~ 1L,
      sex=="Female" ~ 2L,
      TRUE ~ -9L),
    ages,
    agesq,
    centre=assessment_centre,
    bmi,
    neversmoker,
    townsend
  ) %>%
  semi_join(fam, by = c("FID", "IID"))


stopifnot(!any(is.na(covar_plink$FID)))
stopifnot(!any(is.na(covar_plink$IID)))



write.table(
  covar_plink,
  paste0("Outputs/plink_covar.txt"),
  quote = FALSE,
  row.names = FALSE
)



####################

# phenotype file

library(dplyr)

pheno<- read.csv("./Outputs/CVD_Dec_2025.csv") # <- phenotype data frame
name<- "CVD2025"

pheno_names<- names(pheno)[!names(pheno)=="eid"]

pheno_plink <- pheno %>%
  mutate(
    FID = as.numeric(eid),
    IID = as.numeric(eid),
    across(
      all_of(pheno_names),
      ~ if_else(.x, 2L, 1L, missing = -9L)
      )
    ) %>% 
  select(FID, IID, all_of(pheno_names)) %>%
  semi_join(fam, by = c("FID", "IID")) %>%
  semi_join(covar_plink,by = c("FID", "IID"))


stopifnot(!any(is.na(pheno_plink$FID)))
stopifnot(!any(is.na(pheno_plink$IID)))


write.table(
  pheno_plink,
  paste0("Outputs/plink_phenotypes","_",name,".txt"),
  quote = FALSE,
  row.names = FALSE
)


