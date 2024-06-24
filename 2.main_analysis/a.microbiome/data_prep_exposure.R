rm(list = ls())

library(readr)
library(dplyr)
library(readxl)
library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)

exp_data <- as.data.frame(read_xlsx("2a.GWASs_SNPs_EU/Pool de GWAS microbioma EU.xlsx"))
head(exp_data)
str(exp_data)
names(exp_data)

#I firstly obtain the standard error for every SNP as one of the studies did not include it in the summary statistics
exp_data$se <- get_se(exp_data$beta,exp_data$pval)

#Now I calculate the F-statistic of every SNP, if one was lower than 10 it would be excluded
exp_data$fstat <- (abs(exp_data$beta))^2/exp_data$se^2
cat(c(length(which(exp_data$fstat <= 10)),"associations present a F-statistic value of 10 or lower"))
cat(c("The mean value of the F-statistic is",mean(exp_data$fstat)))
#No SNP has an F-stat < 10 and the mean value for the dataset is 33.127

#Checking for duplicates
length(which(duplicated(exp_data$SNP)))
exp_data[which(duplicated(exp_data$SNP)),"SNP"]

exp_data <- exp_data %>%
  group_by(SNP) %>%
  arrange(pval) %>%
  filter(!is.na(SNP)) %>%
  slice(1) %>%
  ungroup()

#Deleting traits with HB as these are for presence and not abundance
exp_data <- exp_data[!grepl("HB",exp_data$Phenotype),]

#Formating the data
exp_data <- format_data(exp_data, type = "exposure")

head(exp_data)

#There is some sort of invisible character in the position column, so I will only keep numbers
exp_data$pos.exposure <- gsub("[^0-9]", "", exp_data$pos.exposure)

write.csv(exp_data,"format_exp_EUR.csv",row.names = FALSE)


#Clumping
exp_data <- read.csv("format_exp_EUR.csv", sep = ";") #Run this if the object has not already been created
str(exp_data)

#Fixing the pval column not being numeric and having problems to turn it because some rows had a comma for decimals instead of a dot
exp_data <- exp_data %>%
  mutate(pval.exposure = gsub(",",".",exp_data$pval.exposure))
exp_data$pval.exposure <- as.numeric(exp_data$pval.exposure)

#Servers Bristol
exp_clump <- ld_clump(
  dplyr::tibble(rsid=exp_data$SNP, 
                pval=exp_data$pval.exposure, 
                id=exp_data$id.exposure),
  clump_kb = 10000, #Clumping kb window
  clump_r2 = 0.001, #if it is >, the SNP with the lowest P-value will be retained
  pop = "EUR")


cat(c("There are",dim(exp_clump)[1],"associations after clumping using rsids and 1000genome ref for EUR"))
write.csv(exp_clump,"exp_clump_EUR.csv", row.names=F)

#Local
exp_clump <- ld_clump(
  dplyr::tibble(rsid=exp_data$SNP, 
                pval=exp_data$pval.exposure, 
                id=exp_data$id),
  clump_kb = 10000, #Clumping kb window
  clump_r2 = 0.001, #if it is >, the SNP with the lowest P-value will be retained
  plink_bin = get_plink_exe(),
  bfile = "C:/bfiles/EUR"
  )

cat(c("There are",dim(exp_clump)[1],"associations after clumping using rsids and 1000genome ref for EUR"))

write.csv(exp_clump,"exp_clump_EUR.csv", row.names=F)

### EAST ASIANS###

exp_data_EAS <- read.csv2("2b.GWAS_SNPs_EA/Microbiome GWAS EA Associations.csv")
head(exp_data_EAS)
names(exp_data_EAS)
str(exp_data_EAS)

exp_data_EAS <- exp_data_EAS %>%
  rename(c("chr" = "Chr",
           "pos" = "Positiona",
           "id" = "Bacterial.group",
           "maf" = "MAFb",
           "se" = "SE",
           "pval" = "P.value"))

#Correcting columns class errors
exp_data_EAS$beta <- as.numeric(exp_data_EAS$beta)
exp_data_EAS$se <- as.numeric(exp_data_EAS$se)
exp_data_EAS$pval <- as.numeric(exp_data_EAS$pval)

str(exp_data_EAS)

#Now I calculate the F-statistic of every SNP, if one was lower than 10 it would be excluded
exp_data_EAS$fstat <- (abs(exp_data_EAS$beta))^2/exp_data_EAS$se^2
cat(c(length(which(exp_data_EAS$fstat <= 10)),"associations present a F-statistic value of 10 or lower"))
cat(c("The mean value of the F-statistic is",mean(exp_data_EAS$fstat)))

#Formatting the data
exp_data_EAS <- format_data(exp_data_EAS, type = "exposure",
                            phenotype_col ="id")
head(exp_data_EAS)

#Deleting every semicolon from the dataset as these mess up the csv file
exp_data_EAS <- exp_data_EAS %>%
  mutate_all(~ gsub(";", "", .))

write.csv(exp_data_EAS, "2b.GWAS_SNPs_EA/format_exp_EAS.csv", row.names = FALSE)

#Check on your dataset (how many assocations before clumping)
dim(exp_data_EAS)[1]
# [1] 84 associations
length(which(!is.na(exp_data_EAS$id.exposure)))*100/length(exp_data_EAS$SNP)
# [1] 100% of the associations with a rsid
length(which(is.na(exp_data_EAS$id.exposure)))*100/length(exp_data_EAS$SNP)
# [1] 0% of the associations with no rsid


#Clumping

#Bristol servers
exp_EAS_clump <- ld_clump(
  dplyr::tibble(rsid=exp_data_EAS$SNP, 
                pval=exp_data_EAS$pval.exposure, 
                id=exp_data_EAS$exposure),
  clump_kb = 10000, #Clumping kb window
  clump_r2 = 0.001, #if it is >, the SNP with the lowest P-value will be retained
  pop = "EAS")

cat(c("There are",dim(exp_EAS_clump)[1],"associations after clumping using rsids and 1000genome ref for EAS"))

write.csv2(exp_EAS_clump,"../2.Selection_SNPs_exposure/3.exp_clump_EAS.csv", row.names=F)

#Local
exp_EAS_clump2 <- ld_clump(
  dplyr::tibble(rsid=exp_data_EAS$id.exposure, 
                pval=exp_data_EAS$pval, 
                id=exp_data_EAS$trait),
  clump_kb = 10000, #Clumping kb window
  clump_r2 = 0.001, #if it is >0.01, the SNP with the lowest P-value will be retained as they will be in LD (if it is <0.01, they both will be retained because they are independent)
  pop = "EAS",
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "C:/Users/alvar/Downloads/1kg.v3.tgz/EAS")

table(exp_EAS_clump$id)
