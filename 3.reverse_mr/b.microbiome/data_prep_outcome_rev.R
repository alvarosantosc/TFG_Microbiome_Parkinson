rm(list = ls())

library(readr)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(readxl)

out_data <- as.data.frame(read_xlsx("2a.GWASs_SNPs_EU/Pool de GWAS microbioma EU.xlsx")) #Pool of selected GWASs summary statistics filtered
head(out_data)
str(out_data)
names(out_data)

#I firstly obtain the standard error for every SNP as one of the studies did not include it in the summary statistics
out_data$se <- get_se(out_data$beta,out_data$pval)

#Formating the data
out_data <- format_data(out_data, type = "outcome")

head(out_data)

#Saving the data
write.csv(out_data,"../3.Selection_SNPs_outcome/rev_MR/format_out_EUR.csv", row.names = FALSE)


### read_outcome_data() ###
exp_data <- read.csv("rev_mr/exp_pool.csv")
exp_clump <- read.csv("rev_mr/exp_clump_pool.csv")
exp_data <- exp_data[which(exp_data$SNP %in% exp_clump$rsid),]

outcome_data <- read_outcome_data(
  snps = exp_data$SNP,
  filename = "../3.Selection_SNPs_outcome/rev_MR/format_out_EUR.csv",
  sep = ";",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  eaf_col = "eaf.outcome",
  pval_col = "pval.outcome",
  id_col = "id.outcome")

#Error in `[[<-.data.frame`(`*tmp*`, type, value = "outcome") : 
#  replacement has 1 row, data has 0

any(out_data$SNP %in% exp_data$SNP) #No results, so the analysis cannot be performed
