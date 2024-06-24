rm(list =ls())


library(readr)
library(dplyr)

IL17_pw <- read_tsv("IL-17 pathway.tsv") #Proteins in the IL17 pathway downloaded from Reactome
IL17_pw <- as.data.frame(IL17_pw)
head(IL17_pw)
names(IL17_pw)
dim(IL17_pw)

#Keeping only the proteins
IL17_pw <- subset(IL17_pw, grepl("Proteins",MoleculeType))
#Deleting a fourth column that shouldn't be there
IL17_pw <- IL17_pw[1:3]

GWAS_proteome <- read_csv2("gwas_proteome.csv", skip = 4) #Summary statistics downloaded from the UK Biobank proteome GWAS

#Keeping only rows with matching identifier with IL-17 pathway
GWAS_IL17 <- GWAS_proteome[GWAS_proteome$`Target UniProt` %in% IL17_pw$Identifier,]
table(GWAS_IL17$`Target UniProt`,useNA = "always")
length(levels(as.factor(GWAS_IL17$`Target UniProt`)))

#35 pQTLs for 13 different proteins

#Getting only cis pQTLs
GWAS_IL17_cis <- as.data.frame(GWAS_IL17[which(GWAS_IL17$`cis/trans` == "cis"),])

#Checking we only have cis SNPs
GWAS_IL17_cis[,c("Target UniProt","rsID","cis/trans")]

#Checking if there are any duplicated SNPs or proteins
length(which(duplicated(GWAS_IL17_cis$rsID)))
length(which(duplicated(GWAS_IL17_cis$Target.UniProt)))
#We have 8 different SNPs, each one for a different protein


write.csv2(GWAS_IL17_cis, "IL-17 Pathway GWAS cis.csv")
