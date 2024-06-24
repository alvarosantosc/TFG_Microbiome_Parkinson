rm(list = ls())


library(readr)
library(dplyr)
library(readxl)
library(openxlsx)

GWAS_1 <- as.data.frame(read_xlsx("GWAS 1 Parkinson_raw.xlsx")) #Summary statistics from selected GWAS 1
head(GWAS_1)
names(GWAS_1)

#As P values are given in a format that R might have problems with, this function will solve it
convert_to_scientific <- function(x) {
  x <- gsub(" × ", "e", x)  # Replace " × " with "e"
  x <- gsub("10−", "-", x) # Replace "10−" with "-"
  return(as.numeric(x))
}

GWAS_1$P.Discovery <- sapply(GWAS_1$P.Discovery, convert_to_scientific)
GWAS_1$P.Replication <- sapply(GWAS_1$P.Replication, convert_to_scientific)
GWAS_1$P.Joint <- sapply(GWAS_1$P.Joint, convert_to_scientific)

#Suponiendo que tengo que utilizar el P.Joint y no el de discovery
GWAS_1_p <- GWAS_1[which(GWAS_1$`P.Joint` < 5e-8),]

length(which(duplicated(GWAS_1_p$SNP)))
#No duplicates

write.xlsx(GWAS_1_p,"GWAS 1 Parkinson filtrado.xlsx")


GWAS_2 <- as.data.frame(read_xlsx("GWAS 2 Parkinson_raw.xlsx")) #Summary statistics from selected GWAS 2
head(GWAS_2)
names(GWAS_2)

GWAS_2$`P joint` <- sapply(GWAS_2$`P joint`, convert_to_scientific)

GWAS_2$`P discovery` <- sapply(GWAS_2$`P discovery`, convert_to_scientific)

GWAS_2$`P NeuroX` <- sapply(GWAS_2$`P NeuroX`, convert_to_scientific)

GWAS_2_p <- GWAS_2[which(GWAS_2$`P joint` < 5e-8),]

length(which(duplicated(GWAS_2_p$SNP)))
#No duplicates

write.xlsx(GWAS_2_p,"GWAS 2 Parkinson filtrado.xlsx")


GWAS_4 <- as.data.frame(read.csv2("GWAS 4 Parkinson_raw.csv")) #Summary statistics from selected GWAS 4
head(GWAS_4)
names(GWAS_4)

GWAS_4$pvalue <- as.numeric(GWAS_4$pvalue)
GWAS_4$OR <- as.numeric(GWAS_4$OR)

GWAS_4_p <- GWAS_4[which(GWAS_4$pvalue < 5e-8),]

length(which(duplicated(GWAS_4_p$assay.name)))
#No duplicates

write.xlsx(GWAS_4_p,"GWAS 4 Parkinson filtrado.xlsx")
