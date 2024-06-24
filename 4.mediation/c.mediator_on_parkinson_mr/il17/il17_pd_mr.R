rm(list = ls())

library(readr)
library(dplyr)
library(TwoSampleMR)
library(ggforestplot)

exp_data <- read.csv("format_IL17_exp.csv", sep = ";")
exp_clump <- read.csv("IL17_clump.csv", sep = ";")
exp_data <- exp_data[which(exp_clump$rsid %in% exp_data$SNP),]

out_data <- read.csv("../../3.Selection_SNPs_outcome/2a.GWASs_SNPs_EU/outcome_data_EUR.csv", sep = ";")

any(out_data$SNP %in% exp_data$SNP)
#None of the 7 pQTLs from the IL-17 pathway are present in the outcome dataset

out_data2 <- read_outcome_data(
  snps = exp_data$SNP,
  filename = "../../3.Selection_SNPs_outcome/2a.GWASs_SNPs_EU/outcome_data_EUR.csv",
  sep = ";",
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

###Extracting outcome data from the package
ao <- available_outcomes()
nalls <- ao[grep("Nalls", ao$author, ignore.case = T),]

outcome <- extract_outcome_data(exp_data$SNP,"ieu-b-7")

harm <- harmonise_data(exposure_dat = exp_data, outcome_dat = outcome, action = 2)

write.csv(harm,"../../4.Harmonize/IL17_Nalls_harm.csv", row.names = FALSE)


#Analysis
harm <- read.csv("../../4.Harmonize/IL17_Nalls_harm.csv") #Run this if object not already created

str(harm)
harm <- harm %>%
  mutate(pval.exposure = gsub(",",".", harm$pval.exposure))
harm$pval.exposure <- as.numeric(harm$pval.exposure)

mr_results <- mr(harm,method_list = c("mr_wald_ratio","mr_ivw"))

mr_results <- generate_odds_ratios(mr_results)

mr_results <- subset(mr_results, select = -c(lo_ci,up_ci))

write.csv2(mr_results, "../../../results/MR/mediation_IL17_PD_Nalls.csv", row.names = FALSE)

###Plotting###
ggforestplot::forestplot(
  df = mr_results,
  name = exposure,
  estimate = b,
  se = se,
  pvalue = pval,
  psignif = 0.05,
  title = "IL17 related proteins associations with Parkinson's Disease risk",
  xlab = "Odds ratio",
  logodds = TRUE
)
