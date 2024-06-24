rm(list = ls())

library(readr)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)

data_GWAS <- read.csv("Icelandic_SNPs_B12.csv", sep = ";")
head(data_GWAS)
names(data_GWAS)
str(data_GWAS)

#Keeping genome-wide associated SNPs, checking for duplicates and obtaining the SE
data_GWAS <- data_GWAS[which(data_GWAS$pval <= 5e-8),]
length(which(duplicated(data_GWAS$SNP)))
data_GWAS$se <- get_se(data_GWAS$beta, data_GWAS$pval)

#Formatting the data
data_GWAS <- mutate(data_GWAS,"Phenotype" = "B12_serum")
data_GWAS <- format_data(data_GWAS, type = "exposure")

###Obtaining studies from ieu##
ao <- available_outcomes()
data.frame(ao[c(grep("b12", ao$trait)),])
data.frame(ao[c(grep("B12", ao$trait)),])

#All of these studies are way lower in sample size to our threshold or their trait is not "Vitamin B12"

data_ieu <- extract_instruments(outcomes = 'ukb-a-135')

write.csv(exp_data, "format_B12_exp.csv", row.names = FALSE)


#Clumping
genetics.binaRies::get_plink_binary()

exp_clump <- ld_clump(
  dplyr::tibble(rsid=exp_data$SNP, 
                pval=exp_data$pval.exposure, 
                id=exp_data$id.exposure),
  clump_kb = 10000, #Clumping kb window
  clump_r2 = 0.001, #if it is >, the SNP with the lowest P-value will be retained
  pop = "EUR")

write.csv(exp_clump,"B12_clump.csv", row.names = FALSE)

exp_data <- exp_data[which(exp_clump$rsid %in% exp_data$SNP)]


###OUTCOME##
out_data <- format_data(data, type = "outcome")

write.csv(out_data,"format_B12_out.csv", row.names = FALSE)
