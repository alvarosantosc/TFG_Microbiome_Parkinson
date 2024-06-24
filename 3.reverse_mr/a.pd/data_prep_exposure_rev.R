rm(list = ls())


library(readr)
library(readxl)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(vroom)
library(plinkbinr)


#Getting the GWAS pool data
exp_data_pool <- as.data.frame(read_xlsx("2a.GWASs_SNPs_EU/Pool de SNPs Parkinson.xlsx")) #Summary statistics from selected GWASs from GWAS catalog
head(exp_data_pool)
names(exp_data_pool)
str(exp_data_pool)

exp_data_pool <- exp_data_pool %>%
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
head(exp_data_pool)

#Calculating the missing values for standard error and beta value
exp_data_pool$or <- as.numeric(exp_data_pool$or)
exp_data_pool$beta <- log(exp_data_pool$or)
exp_data_pool$se <- get_se(exp_data_pool$beta,exp_data_pool$pval)

exp_data_pool$eaf <- as.numeric(exp_data_pool$eaf)

#Looking for duplicates to keep only one
exp_data_pool[which(duplicated(exp_data_pool$SNP)),c("SNP")]
exp_data_pool[which(exp_data_pool$SNP=="rs356182"),]
exp_data_pool[which(exp_data_pool$SNP=="rs11158026"),]

exp_data_pool <- exp_data_pool %>%
  group_by(SNP) %>%
  arrange(pval) %>%
  slice(1) %>%
  ungroup()

exp_data_pool$Phenotype <- "Parkinson's Disease"
exp_data_pool$id <- "GWAS Catalog"
exp_data_pool <- format_data(exp_data_pool, type = "exposure")
###Getting the Nalls 2019 dataset###
ao <- available_outcomes()

#a <- tophits(id="ieu-b-7", clump=1) #Devuelve los SNPs con un pval <= 5e-8 y eliges si hace clumping o no
exp_data_nalls <- extract_instruments(outcomes = "ieu-b-7", p1 = 1) #Hace lo mismo que tophits pero por defecto hace clumping y da los datos formateados para exposure, así que mejor esta función

exp_data_nalls <- subset(exp_data_nalls, select = -c(samplesize.exposure, data_source.exposure))

exp_data_nalls$exposure <- "Parkinson's Disease" #So that the clumping step understands that it is always the same exposure
###Getting the FinnGen dataset###
exp_data_finngen <- data.frame(vroom("../1.Selection_GWAS/Parkinson/summary_stats_finngen_R10_G6_PARKINSON.gz"))

exp_data_finngen <- exp_data_finngen[which(exp_data_finngen$pval <= 5e-8),]

cat(c("There are",length(which(duplicated(exp_data_finngen$rsids))),"duplicated rsIDs"))

exp_data_finngen[which(duplicated(exp_data_finngen$rsids)),]
exp_data_finngen[which(exp_data_finngen$rsids == "rs2696531"),]
exp_data_finngen[which(exp_data_finngen$rsids == "rs2696530"),]

#Filtering by p-val only if this differs by more than 1 magnitude order
exp_data_finngen <- as.data.frame(exp_data_finngen %>%
                                   group_by(rsids) %>%
                                   mutate(highest_p_value = min(pval)) %>%
                                   filter(pval / highest_p_value < 10) %>%
                                   select(-highest_p_value) %>%
                                   ungroup())

cat(c(dim(exp_data_finngen)[1] - dim(exp_data_finngen)[1],"duplicates have been deleted"))
#Filtering by beta value
exp_data_finngen <- as.data.frame(exp_data_finngen %>%
                                   group_by(rsids) %>%
                                   arrange(desc(abs(beta))) %>%
                                   slice(1) %>%
                                   ungroup())

exp_data_finngen <- exp_data_finngen %>%
  rename(c("chr" = "X.chrom",
           "other_allele" = "ref",
           "effect_allele" = "alt",
           "SNP" = "rsids",
           "se" = "sebeta",
           "eaf" = "af_alt"))

#I have identified several rsIDs that are not SNPs but indels, keeping only SNPs
exp_data_finngen <- exp_data_finngen %>%
  filter(nchar(other_allele) == 1) %>%
  filter(nchar(effect_allele) == 1)

#Adding the phenotype and formating
exp_data_finngen$Phenotype <- "Parkinson's Disease"
exp_data_finngen$id <- "Finngen"
exp_data_finngen <- format_data(exp_data_finngen, type = "exposure")



#Unifying the three datasets
exp_data <- rbind(exp_data_nalls,exp_data_pool,exp_data_finngen)

#Checking duplicates among sources
length(which(duplicated(exp_data$SNP)))
exp_data[which(duplicated(exp_data$SNP)),"SNP"]
exp_data[which(exp_data$SNP == "rs6430538"),c("SNP","pval.exposure","beta.exposure","id.exposure")]
exp_data[which(exp_data$SNP == "rs17649553"),c("SNP","pval.exposure","beta.exposure","id.exposure")]
exp_data[which(exp_data$SNP == "rs2230288"),c("SNP","pval.exposure","beta.exposure","id.exposure")]
exp_data[which(exp_data$SNP == "rs73984689"),c("SNP","pval.exposure","beta.exposure","id.exposure")]
exp_data[which(exp_data$SNP == "rs365825"),c("SNP","pval.exposure","beta.exposure","id.exposure")]

#They all differ widely in the pval, so only this parameter and not beta will be considered
exp_data <- as.data.frame(exp_data %>%
  group_by(SNP) %>%
  arrange((abs(pval.exposure))) %>%
  slice(1) %>%
  ungroup())

write.csv(exp_data,"../2.Selection_SNPs_exposure/rev_mr/exp_pool.csv", row.names = FALSE)

###Clumping###
exp_clump <- ld_clump(
  dplyr::tibble(rsid=exp_data$SNP, pval=exp_data$pval.exposure, id=exp_data$exposure),
  clump_kb = 10000, #Clumping kb window
  clump_r2 = 0.001, #if it is >, the SNP with the lowest P-value will be retained
  pop = "EUR")

#Local
exp_clump <- ld_clump(
  dplyr::tibble(rsid=exp_data$SNP, 
                pval=exp_data$pval.exposure, 
                id=exp_data$exposure),
  clump_kb = 10000, #Clumping kb window
  clump_r2 = 0.001, #if it is >, the SNP with the lowest P-value will be retained
  plink_bin = get_plink_exe(),
  bfile = "C:/bfiles/EUR"
)

cat(c("There are",dim(exp_clump)[1],"associations after clumping using rsids and 1000genome ref for EUR"))

write.csv(exp_clump,"../2.Selection_SNPs_exposure/rev_mr/exp_clump_pool.csv", row.names = FALSE)

