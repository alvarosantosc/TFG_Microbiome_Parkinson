rm(list = ls())

library(readr)
library(dplyr)
library(TwoSampleMR)

exp_data <- read.csv("../../2.Selection_SNPs_exposure/2a.GWASs_SNPs_EU/format_exp_EUR.csv")
exp_clump <- read.csv("../../2.Selection_SNPs_exposure/exp_clump_EUR.csv")
exp_data <- exp_data[which(exp_clump$rsid %in% exp_data$SNP),]


out_data <- read.csv("format_B12_out.csv")

any(out_data$SNP %in% exp_data$SNP)
#None of the 8 SNPs from the B12 GWAS are present in the microbiome dataset

out_data2 <- read_outcome_data(
  snps = exp_data$SNP,
  filename = "format_B12_out.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  eaf_col = "eaf.outcome",
  pval_col = "pval.outcome",
  phenotype_col = "outcome"
)

#Error in `$<-.data.frame`(`*tmp*`, "mr_keep.outcome", value = TRUE) : 
#  replacement has 1 row, data has 0
