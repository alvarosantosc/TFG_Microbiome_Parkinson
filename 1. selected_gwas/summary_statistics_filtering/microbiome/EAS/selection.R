rm(list = ls())


library(readr)
library(dplyr)
library(tidyr)

GWAS_associations <- as.data.frame(read_csv2("GWAS Microbiome EA associations_raw.csv", skip = 2, n_max = 85)) #Summary statistics from selected GWASs
head(GWAS_associations)
names(GWAS_associations)
dim(GWAS_associations)


#Here I separate into 2 columns the one containing both the effect size and the standard error
GWAS_associations <- GWAS_associations %>%
  separate(`Effect size (SE)c`, into = c("beta", "SE"), sep = " \\(") %>%
  mutate(SE = substr(SE, 1, nchar(SE) - 1)) %>%
  mutate(beta = as.numeric(beta),
         SE = as.numeric(SE))

GWAS_associations <- GWAS_associations %>%
  separate(`Minor/major allele`, into = c("effect_allele","other_allele", sep = "/")) %>%
  select(-`/`)

length(which(duplicated(GWAS_associations$SNP)))
#There is only 1 SNP repeated
#There must be for sure an easier way to print duplicated SNPs and their bacteria, but credits to GPT for this one
GWAS_associations %>%
  filter(duplicated(SNP) | duplicated(SNP, fromLast = TRUE)) %>%
  select(SNP, `Bacterial group`, beta) %>%
  print()
#So the only duplicated SNP comes from different bacteria

#Filtering for the SNP with the highest absolute beta value
GWAS_associations_filter <- as.data.frame(GWAS_associations %>%
  group_by(SNP) %>%
  arrange(desc(abs(beta))) %>%
  slice(1) %>%
  ungroup())


length(which(duplicated(GWAS_associations_filter$SNP)))
#No duplicated SNPs anymore
GWAS_associations_filter[which(GWAS_associations_filter$SNP == "rs1448221"), c("SNP","beta")]
#The remaining SNP is the one with the highest absolute beta value

write.csv(GWAS_associations_filter,"Microbiome GWAS EA Associations.csv", row.names = FALSE)
