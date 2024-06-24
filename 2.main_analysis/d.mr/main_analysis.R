rm(list=ls())

library(readr)
library(TwoSampleMR)


#ANALYSIS#
harm <- read.csv("Nalls_harm_main.csv")
str(harm)

harm <- harm %>%
  mutate(pval.exposure = gsub(",",".", harm$pval.exposure))
harm$pval.exposure <- as.numeric(harm$pval.exposure)


mr_results <- mr(harm,method_list = c("mr_wald_ratio","mr_ivw"))

mr_results <- generate_odds_ratios(mr_results)

mr_results <- subset(mr_results, select = -c(lo_ci,up_ci))

#Saving the results
write.csv2(mr_results, "../../results/MR/main_analysis_Nalls.csv", row.names = FALSE)

mr_results <- read.csv2("../../results/MR/main_analysis_Nalls.csv") #Run this is the object is not already created

mr_results_sign <- mr_results[which(mr_results$pval <= 0.05),]
#mr_results_sign <- mr_results_sign[which(mr_results_sign$or < 10),]

write.csv2(mr_results_sign, "../../results/MR/main_analysis_Nalls_significative.csv", row.names = FALSE)

