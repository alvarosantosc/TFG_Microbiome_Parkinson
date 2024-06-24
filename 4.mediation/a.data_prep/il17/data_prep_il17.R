rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(TwoSampleMR)
library(ieugwasr)

data <- read.csv2("IL-17 Pathway GWAS cis.csv")
head(data)
names(data)
str(data)

#Extracting the effect allele (A1) and other allele (A0)
data <- data %>% 
  separate(Variant.ID..CHROM.GENPOS..hg37..A0.A1.imp.v1., into = c("CHROM", "GENPOS_hg37", "A0", "A1", "imp", "v1"), sep = ":")

data <- data %>%
  rename(c("chr" = "CHROM",
           "pos" = "GENPOS..hg38.",
           "other_allele" = "A0",
           "effect_allele" = "A1",
           "eaf" = "A1FREQ..discovery.",
           "Phenotype" = "Assay.Target",
           "id" = "Target.UniProt",
           "SNP" = "rsID",
           "beta" = "BETA..replication.",
           "se" = "SE..replication.",
           "pval" = "log10.p...replication."))

#One of the pQTLs is not a SNP but an indel
data <- data %>%
  filter(nchar(other_allele) == 1) %>%
  filter(nchar(effect_allele) == 1)

#Transforming pval from log scale
data$pval <- 10^(-data$pval)

#Formatting the data
exp_data <- format_data(data, type = "exposure")

write.csv(exp_data, "format_IL17_exp.csv", row.names = FALSE)


#Clumping
genetics.binaRies::get_plink_binary()

exp_clump <- ld_clump(
  dplyr::tibble(rsid=exp_data$SNP, 
                pval=exp_data$pval.exposure, 
                id=exp_data$id.exposure),
  clump_kb = 10000, #Clumping kb window
  clump_r2 = 0.001, #if it is >, the SNP with the lowest P-value will be retained
  pop = "EUR")

write.csv(exp_clump, "IL17_clump.csv", row.names = FALSE)



###OUTCOME##
out_data <- format_data(data, type = "outcome")

write.csv(out_data,"format_IL17_out.csv", row.names = FALSE)
