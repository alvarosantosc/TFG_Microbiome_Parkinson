#########################
#R environment
#########################

rm(list = ls())

#Load libraries
library(readr)
library(dplyr)



#########################
#Read data
#########################

GWAS_B12 <- read_tsv("GWAS Vitamin B12_raw.tsv") #Studies downloaded from GWAS Catalog with EFO_0004620
GWAS_B12<- as.data.frame(GWAS_B12)
head(GWAS_B12)
names(GWAS_B12)
dim(GWAS_B12)



#########################
#QC 
#########################

#Check if there are duplicated studies

# GWAS_B12 <- GWAS_B12 %>%
#   distinct(.keep_all = TRUE)
# dim(GWAS_B12)

GWAS <- GWAS_B12 %>%
  filter(reportedTrait == "Vitamin B12 levels") 
#AFS: you are excluding "Cobalamin (Vitamin B12) [Mass/volume] in Serum or Plasma"
dim(GWAS_B12)



#########################
#Selection of EU GWAS
#########################

GWAS_EU <- subset(GWAS, grepl("European",discoverySampleAncestry))
dim(GWAS_EU)

#No need to filter more as we already have only europeans with n>1000

#Add column with sample size
GWAS_EU$sample_size <- as.numeric(gsub("\\D", "", GWAS_EU$discoverySampleAncestry))
table(GWAS_EU$sample_size)

#Order by sample size
head(GWAS_EU)
tail(GWAS_EU)
GWAS_EU <- GWAS_EU[order(GWAS_EU$sample_size, decreasing = T),]
head(GWAS_EU)
tail(GWAS_EU)

#Check on potential duplicated studies 

#Save dataset
GWAS_EU <- GWAS_EU[,c("firstAuthor","publicationDate","journal","title","pubmedId",                  
                      "reportedTrait",
                      "discoverySampleAncestry","initialSampleDescription","sample_size",
                      "replicationSampleAncestry","replicateSampleDescription")]
write.table(GWAS_EU,"B12_GWAS.csv", sep = "\t", row.names = F) 
# 
# #Keeping only studies with n>1000
# GWAS_EU_minN100000 <- as.data.frame(subset(GWAS_EU, sample_size >= 100000))
# table(GWAS_EU_minN1000$sample_size)
# 
# write.table(GWAS_EU_minN1000,"B12_GWAS_N1000.csv", sep = "\t", row.names = F) 

#Save IDs of the selected GWAS
GWAS_EU_ids <- GWAS_EU[!duplicated(GWAS_EU$pubmedId),c("firstAuthor","publicationDate","journal","title","pubmedId")]
write.table(GWAS_EU_ids,"B12_GWAS_studies.csv", sep = "\t", row.names = F) 

#Save exposures of the selected GWAS 
length(levels(as.factor(GWAS_EU$reportedTrait)))
length(unique(levels(as.factor(GWAS_EU$reportedTrait))))
EUR_traits <- GWAS_EU[,c("reportedTrait")]
write.table(EUR_traits,"B12_traits.csv", sep = "\t", row.names = F, col.names = F) 

