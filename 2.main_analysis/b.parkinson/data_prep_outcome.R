rm(list = ls())


library(readr)
library(readxl)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(vroom)

#sumstats from Finngen
# setwd("C:/Users/rmjlaf0/OneDrive - University College London/General - ALBA_ALVARO/Project/3.Analysis/data/1.Selection_GWAS/Parkinson")
# parkinson <- vroom("summary_stats_finngen_R10_G6_PARKINSON.gz")


###Unifying PD data from GWAS Catalog, FinnGen and available outcomes###

#Getting the GWAS pool data
out_data_pool <- as.data.frame(read_xlsx("2a.GWASs_SNPs_EU/Pool de SNPs Parkinson.xlsx")) #Pool of every extracted SNP from the GWAS catalog studies
head(out_data_pool)
names(out_data_pool)
str(out_data_pool)

out_data_pool <- out_data_pool %>%
  rename(c("SNP" = "rsID",
           "chr" = "chromosome",
           "pos" = "position",
           "eaf" = "EAF",
           "maf" = "MAF",
           "beta" = "BETA",
           "se" = "SE",
           "pval" = "pVal",
           "or" = "OR",
           "id" = "reportedTrait"
  ))
head(out_data_pool)

#Calculating the missing values for standard error and beta value
out_data_pool$or <- as.numeric(out_data_pool$or)
out_data_pool$beta <- log(out_data_pool$or)
out_data_pool$se <- get_se(out_data_pool$beta,out_data_pool$pval)

out_data_pool$eaf <- as.numeric(out_data_pool$eaf)

#Looking for duplicates to keep only one
out_data_pool[which(duplicated(out_data_pool$SNP)),c("SNP")]
out_data_pool[which(out_data_pool$SNP=="rs356182"),]
out_data_pool[which(out_data_pool$SNP=="rs11158026"),]

out_data_pool <- out_data_pool %>%
  group_by(SNP) %>%
  arrange(pval) %>%
  slice(1) %>%
  ungroup()

out_data_pool$Phenotype <- "Parkinson's Disease"
out_data_pool$id <- "GWAS Catalog"
out_data_pool <- format_data(out_data_pool, type = "outcome")


###Getting the Nalls 2019 dataset###
ao <- available_outcomes()

out_data_nalls <- tophits(id="ieu-b-7", clump = 0) #Aquí usamos tophits para obtener los SNPs porque extract_instruments lo devolvería clumped y formateado para exposure

out_data_nalls <- out_data_nalls %>%
  rename(c("pos" = "position",
           "pval" = "p",
           "effect_allele" = "ea",
           "other_allele" = "nea",
           "SNP" = "rsid"))

out_data_nalls$Phenotype <- "Parkinson's Disease"

#Checking for duplicates
length(which(duplicated(out_data_nalls$SNP)))
out_data_nalls[which(duplicated(out_data_nalls$SNP)),"SNP"]
#Filtering by p-val only if this differs by more than 1 magnitude order
out_data_nalls <- as.data.frame(out_data_nalls %>%
                                    group_by(SNP) %>%
                                    mutate(highest_p_value = min(pval)) %>%
                                    filter(pval / highest_p_value < 10) %>%
                                    select(-highest_p_value) %>%
                                    ungroup())

#Filtering by beta value
out_data_nalls <- as.data.frame(out_data_nalls %>%
                                    group_by(SNP) %>%
                                    arrange(desc(abs(beta))) %>%
                                    slice(1) %>%
                                    ungroup())

out_data_nalls <- format_data(out_data_nalls, type = "outcome")


###Getting the FinnGen dataset###
out_data_finngen <- data.frame(vroom("../1.Selection_GWAS/Parkinson/summary_stats_finngen_R10_G6_PARKINSON.gz"))

out_data_finngen <- out_data_finngen[which(out_data_finngen$pval <= 5e-8),]

cat(c("There are",length(which(duplicated(out_data_finngen$rsids))),"duplicated rsIDs"))

out_data_finngen[which(duplicated(out_data_finngen$rsids)),]
out_data_finngen[which(out_data_finngen$rsids == "rs2696531"),]
out_data_finngen[which(out_data_finngen$rsids == "rs2696530"),]

#Filtering by p-val only if this differs by more than 1 magnitude order
out_data_finngen <- as.data.frame(out_data_finngen %>%
                                    group_by(rsids) %>%
                                    mutate(highest_p_value = min(pval)) %>%
                                    filter(pval / highest_p_value < 10) %>%
                                    select(-highest_p_value) %>%
                                    ungroup())

#Filtering by beta value
out_data_finngen <- as.data.frame(out_data_finngen %>%
                                    group_by(rsids) %>%
                                    arrange(desc(abs(beta))) %>%
                                    slice(1) %>%
                                    ungroup())

out_data_finngen <- out_data_finngen %>%
  rename(c("chr" = "X.chrom",
           "other_allele" = "ref",
           "effect_allele" = "alt",
           "SNP" = "rsids",
           "se" = "sebeta",
           "eaf" = "af_alt"))

#I have identified several rsIDs that are not SNPs but indels, keeping only SNPs
out_data_finngen <- out_data_finngen %>%
  filter(nchar(other_allele) == 1) %>%
  filter(nchar(effect_allele) == 1)

#Adding the phenotype and formating
out_data_finngen$Phenotype <- "Parkinson's Disease"
out_data_finngen$id <- "Finngen"
out_data_finngen <- format_data(out_data_finngen, type = "outcome")


#Unifying the three datasets
out_data <- rbind(out_data_nalls,out_data_pool,out_data_finngen)

#Checking for duplicates among studies
cat("There are",length(which(duplicated(out_data$SNP)))," duplicated SNPs")
head(out_data[which(duplicated(out_data$SNP)),"SNP"])

#Filtering by p-val only if this differs by more than 1 magnitude order
out_data <- as.data.frame(out_data %>%
                                    group_by(SNP) %>%
                                    mutate(highest_p_value = min(pval.outcome)) %>%
                                    filter(pval.outcome / highest_p_value < 10) %>%
                                    select(-highest_p_value) %>%
                                    ungroup())

cat(length(which(duplicated(out_data$SNP))),"duplicated SNPs remain")

#Filtering by beta value
out_data <- as.data.frame(out_data %>%
                                    group_by(SNP) %>%
                                    arrange(desc(abs(beta.outcome))) %>%
                                    slice(1) %>%
                                    ungroup())

cat(length(which(duplicated(out_data$SNP))),"duplicated SNPs remain")


write.csv(out_data,"2a.GWASs_SNPs_EU/outcome_data_EUR.csv", row.names = FALSE)




### EAST ASIAN ###

out_data_EAS <- as.data.frame(read.csv2("2b.GWASs_SNPs_EA/GWAS 1 Parkinson EA filtrado.csv"))
head(out_data_EAS)
names(out_data_EAS)
str(out_data_EAS)


#Renaming columns
out_data_EAS <- out_data_EAS %>%
  rename(c("chr" = "CHR",
           "pos" = "BP..11219.",
           "left_gene" = "Left.gene",
           "right_gene" = "Right.gene",
           "effect_allele" = "A1",
           "pval" = "Meta.P",
           "or" = "Meta..OR.discovery...validation1",
           "outcome" = "reportedTrait"))

#Fixing columns class errors
str(out_data_EAS)

out_data_EAS$pval <- as.numeric(out_data_EAS$pval)
out_data_EAS$or <- as.numeric(out_data_EAS$or)

#Calculating beta and standard error
out_data_EAS$beta <- log(out_data_EAS$or)
out_data_EAS$se <- get_se(out_data_EAS$beta,out_data_EAS$pval)

#No need to check for duplicates as only 1 study is being used and it was already checked when obtaining its summary statistics

#Formating the data
out_data_EAS <- format_data(out_data_EAS, type = "outcome")


#Saving the data
write.csv(out_data_EAS,"2b.GWASs_SNPs_EA/outcome_data_EAS.csv", row.names = FALSE)

### Extract_outcome_data() ###

exp_clump_EAS <- read.csv("C:/Users/alvar/OneDrive - University College London/General - ALBA_ALVARO/Project/3.Analysis/data/2.Selection_SNPs_exposure/3.exp_clump_EAS.csv", sep = ";")
exp_data_EAS <- read.csv("../2.Selection_SNPs_exposure/2b.GWAS_SNPs_EA/format_exp_EAS.csv", sep = ";")
exp_data_EAS <- exp_data_EAS[which(exp_data_EAS$SNP %in% exp_clump_EAS$rsid),]

outcome_EAS <- read_outcome_data(
  snps = exp_data_EAS$SNP,
  filename = "2b.GWASs_SNPs_EA/outcome_data_EAS.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  pval_col = "pval.outcome",
  phenotype_col = "outcome")
