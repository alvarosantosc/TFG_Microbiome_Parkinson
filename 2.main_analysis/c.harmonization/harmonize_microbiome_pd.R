rm(list=ls())

library(TwoSampleMR)
library(readr)

#exposure (clumped)
exp_data <- data.frame(read.csv("2.Selection_SNPs_exposure/format_exp_EUR.csv", sep = ";"))
exp_clumped <- data.frame(read_delim("2.Selection_SNPs_exposure/exp_clump_EUR.csv")) #the dataset of clumped exposure data needs to have all columns

exp_data <- exp_data[which(exp_clumped$rsid %in% exp_data$SNP),]

###outcome
#From GWAS Catalog + FinnGen + tophits Nalls
outcome <- read_outcome_data(
  snps = exp_data$SNP,
  filename = "3.Selection_SNPs_outcome/2a.GWASs_SNPs_EU/outcome_data_EUR.csv",
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
#No outcome data is obtained by this method

#From IEUopenGWAS
ao <-available_outcomes()
nalls <- ao[grep("Nalls", ao$author, ignore.case = T),]

outcome <- extract_outcome_data(exp_data$SNP,"ieu-b-7")

harm <- harmonise_data(exposure_dat = exp_data, outcome_dat = outcome, action = 2)

write.csv(harm,"4.Harmonize/Nalls_harm_main.csv", row.names = FALSE)


#Checking these results
min(outcome$pval.outcome)

#In fact, the lowest pval is just 0.000312997, far from considered in the genome-wide association threshold

