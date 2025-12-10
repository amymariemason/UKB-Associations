####################################################
# create instrument 
####################################################


#  new alcohol instrument define

alcovars99_updated <- read_csv("Inputs/alcovars99_updated.csv")

library(readxl)
Saunders_2022_Supplement <- read_excel("Inputs/Saunders_2022_Supplement.xlsx", 
                                         +     sheet = "alcohol_smoking")
library(tidyverse)

####
old_vars <- alcovars99_updated   %>% filter (Chr!="") %>%
  mutate(chr= Chr, pos=hg38_pos, phenotype="old alcohol") %>%
 select(chr, pos, snp, phenotype)

new_vars<- Saunders_2022_Supplement %>%
  mutate(chr= as.numeric(sub("chr", "", chr)), pos=position_b38) %>% 
  select(chr, pos, snp, phenotype)

both_sets<- rbind(old_vars, new_vars)
dups <- duplicated(both_sets %>% select(chr, pos))
both_sets_unique<- both_sets[!dups,] %>% arrange(chr, pos)

###############################

# use ukbrapr to make instrument

# install current version
library(pacman)
p_load_gh("lcpilling/ukbrapR")


varlist <- data.frame(rsid=both_sets_unique$snp, chr=both_sets_unique$chr) 

imputed_genotypes <- extract_variants(varlist,out_bed = "alcohol_and_smoking")
#> ~1h 40 minutes

dim(imputed_genotypes)


