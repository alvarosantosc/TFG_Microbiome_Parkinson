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

GWAS_PD <- read_tsv("GWAS Parkinson_raw.tsv") #Downloaded studies from GWAS catalog MONDO_0005180
GWAS_PD<- as.data.frame(GWAS_PD)
head(GWAS_PD)
names(GWAS_PD)
dim(GWAS_PD)



#########################
#QC 
#########################


#Select GWAS on the risk of general PD only
GWAS <- GWAS_PD %>%
  filter(reportedTrait == "Parkinson's disease" |reportedTrait == "Parkinsons disease" | reportedTrait == "Parkinson's disease (PheCode 332)" | reportedTrait == "ICD10 G20: parkinsons disease"| reportedTrait == "ICD10 G20: Parkinson's disease") 
###AFS: Add the other PD general (those with IDs, but reflecting risk of general PD)
dim(GWAS)



#########################
#Selection of EU GWAS
#########################

#Subset for EU
GWAS_EU <- subset(GWAS, grepl("European", discoverySampleAncestry))
dim(GWAS_EU)
table(GWAS_EU$pubmedId, useNA = "always")

dim(GWAS_EU)
table(GWAS_EU$pubmedId, useNA = "always")

#Add column with sample size
GWAS_EU$sample_size <- ifelse(
  !is.na(GWAS_EU$replicationSampleAncestry), 
  as.numeric(gsub("^(\\d+) European.*", "\\1", GWAS_EU$discoverySampleAncestry)) + as.numeric(gsub("^(\\d+) European.*", "\\1", GWAS_EU$replicationSampleAncestry)), 
  as.numeric(gsub("^(\\d+) European.*", "\\1", GWAS_EU$discoverySampleAncestry))
)

#By doing this, the sample size of 3 studies is lost, but their sample size is too low to be considered


table(GWAS_EU$sample_size)

#Order by sample size
head(GWAS_EU)
tail(GWAS_EU)
GWAS_EU <- GWAS_EU[order(GWAS_EU$sample_size, decreasing = T),]
head(GWAS_EU)
tail(GWAS_EU)

#Keeping only studies with n>100000
GWAS_EU_minN100000 <- as.data.frame(subset(GWAS_EU, sample_size >= 100000))
table(GWAS_EU_minN100000$sample_size)

#Check on potential duplicated studies 
GWAS_EU_minN100000[which(duplicated(GWAS_EU_minN100000$pubmedId)),c("pubmedId","discoverySampleAncestry")]
GWAS_EU_minN100000 <- GWAS_EU_minN100000 %>%
  group_by(pubmedId) %>%
  arrange(desc(sample_size)) %>%
  slice(1) %>%
  ungroup()

length(which(duplicated(GWAS_EU_minN100000$pubmedId)))
#Save dataset
GWAS_EU_filtered <- GWAS_EU_minN100000[,c("firstAuthor","publicationDate","journal","title","pubmedId",                  
                      "reportedTrait",
                      "discoverySampleAncestry","initialSampleDescription","sample_size",
                      "replicationSampleAncestry","replicateSampleDescription")]
write.table(GWAS_EU_filtered,"Europeans_PD_GWAS.csv", sep = "\t", row.names = F)



#write.table(GWAS_EU_minN1000,"Europeans_PD_GWAS_N1000.csv", sep = "\t", row.names = F) 

#Save IDs of the selected GWAS
GWAS_EU_minN1000_ids <- GWAS_EU_minN1000[!duplicated(GWAS_EU_minN1000$pubmedId),c("firstAuthor","publicationDate","journal","title","pubmedId")]
write.table(GWAS_EU_minN1000_ids,"Europeans_PD_GWAS_studies.csv", sep = "\t", row.names = F) 

#Save exposures of the selected GWAS 
length(levels(as.factor(GWAS_EU_minN100000$reportedTrait)))
length(unique(levels(as.factor(GWAS_EU_minN100000$reportedTrait))))
EUR_traits <- GWAS_EU_minN100000[,c("reportedTrait")]
write.table(EUR_traits,"Europeans_PD_traits.csv", sep = "\t", row.names = F, col.names = F) 

# write.csv2(GWAS_EU_minN1000,"Parkinson GWAS EU.csv", row.names = FALSE)



#########################
#Selection of EA GWAS
#########################

#Subset for EA
GWAS_EA <- subset(GWAS, grepl("East Asian", discoverySampleAncestry))
dim(GWAS_EU)
table(GWAS_EA$pubmedId, useNA = "always")

dim(GWAS_EA)
table(GWAS_EA$pubmedId, useNA = "always")

#Add column with sample size
GWAS_EA$sample_size <- ifelse(
  !is.na(GWAS_EA$replicationSampleAncestry), 
  as.numeric(sub(".*?(\\d+) East Asian.*", "\\1", GWAS_EA$discoverySampleAncestry)) + as.numeric(sub(".*?(\\d+) East Asian.*", "\\1", GWAS_EA$replicationSampleAncestry)), 
  as.numeric(sub(".*?(\\d+) East Asian.*", "\\1", GWAS_EA$discoverySampleAncestry))
)

table(GWAS_EA$sample_size)

#Order by sample size
head(GWAS_EA)
tail(GWAS_EA)
GWAS_EA <- GWAS_EA[order(GWAS_EA$sample_size, decreasing = T),]
head(GWAS_EA)
tail(GWAS_EA)

#Keeping only studies with n>10000
GWAS_EA_minN10000 <- as.data.frame(subset(GWAS_EA, sample_size >= 10000))
table(GWAS_EA_minN10000$sample_size)

#Check on potential duplicated studies 
GWAS_EA_minN10000[which(duplicated(GWAS_EA_minN10000$pubmedId)),c("pubmedId","discoverySampleAncestry")]
GWAS_EA_minN10000 <- GWAS_EA_minN10000 %>%
  group_by(pubmedId) %>%
  arrange(desc(sample_size)) %>%
  slice(1) %>%
  ungroup()

length(which(duplicated(GWAS_EA_minN10000$pubmedId)))

#Save dataset
GWAS_EA_minN10000 <- GWAS_EA_minN10000[,c("firstAuthor","publicationDate","journal","title","pubmedId",                  
                      "reportedTrait",
                      "discoverySampleAncestry","initialSampleDescription","sample_size",
                      "replicationSampleAncestry","replicateSampleDescription")]
write.table(GWAS_EA_minN10000,"EA_PD_GWAS.csv", sep = "\t", row.names = F) 


#Save IDs of the selected GWAS
GWAS_EA_minN10000_ids <- GWAS_EA_minN10000[!duplicated(GWAS_EA_minN10000$pubmedId),c("firstAuthor","publicationDate","journal","title","pubmedId")]
write.table(GWAS_EA_minN10000_ids,"EA_PD_GWAS_studies.csv", sep = "\t", row.names = F) 

#Save exposures of the selected GWAS 
length(levels(as.factor(GWAS_EA_minN10000$reportedTrait)))
length(unique(levels(as.factor(GWAS_EA_minN10000$reportedTrait))))
EA_traits <- GWAS_EA_minN10000[,c("reportedTrait")]
write.table(EA_traits,"EA_PD_traits.csv", sep = "\t", row.names = F, col.names = F) 
