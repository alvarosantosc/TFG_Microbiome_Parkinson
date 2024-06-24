rm(list = ls())

library(readr)
library(dplyr)
library(readxl)

GWAS_1 <- as.data.frame(read_xlsx("GWAS 1 Parkinson_raw.xlsx")) #Summary statistics from selected GWAS
head(GWAS_1)
names(GWAS_1)
str(GWAS_1)

#Fixing columns class erorrs
GWAS_1$`Meta-P` <- as.numeric(GWAS_1$`Meta-P`)
GWAS_1_p <- GWAS_1[which(GWAS_1$`Meta-P` < 5e-8),]

length(which(duplicated(GWAS_1_p$SNP)))

write.csv(GWAS_1_p,"GWAS 1 Parkinson EA filtrado.csv", row.names = FALSE)
