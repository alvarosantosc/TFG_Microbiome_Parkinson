#########################
#R environment
#########################

rm(list = ls())

getwd()
setwd("C:/Users/alvar/OneDrive - University College London/General - ALBA_ALVARO/Project/3.Analysis/data/1.Selection_GWAS/Microbiome/1.EFO_0007874")
# setwd("C:/Users/rmjlaf0/OneDrive - University College London/General - ALBA_ALVARO/Project/3.Analysis/data/1.Selection_GWAS/Microbiome/1.EFO_0007874/")

#Load libraries
library(readr)
library(dplyr)



#########################
#Read data
#########################

GWAS_microbiome <- read_tsv("GWAS microbiome_raw.tsv")
GWAS_microbiome <- as.data.frame(GWAS_microbiome)
head(GWAS_microbiome)
names(GWAS_microbiome)
dim(GWAS_microbiome)



#########################
#QC 
#########################

#Check if there are duplicated studies

###I think this is supposed to eliminate the duplicates, but it did not reduce the number of rows
GWAS <- GWAS_microbiome 
# %>% 
#   distinct(.keep_all = TRUE) #afs: done it below
rm(GWAS_microbiome)
dim(GWAS)

length(levels(as.factor(GWAS$reportedTrait)))
length(GWAS$reportedTrait) - length(levels(as.factor(GWAS$reportedTrait)))
# [1] 10 gut microbiome measures appear more than once
length(which(duplicated(GWAS$reportedTrait)))
table(GWAS$reportedTrait, useNA = "always")
GWAS[which(duplicated(GWAS$reportedTrait)),c("reportedTrait", "pubmedId")]
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[1],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[2],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[3],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[4],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[5],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[6],c("reportedTrait")]), #this trait comes from 3 different studies
     c("reportedTrait", "pubmedId")]
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[7],c("reportedTrait")]), #this trait comes from 3 different studies
     c("reportedTrait", "pubmedId")]
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[8],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[9],c("reportedTrait")]), #this trait comes from 3 different studies
     c("reportedTrait", "pubmedId")]
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[10],c("reportedTrait")]), #this trait comes from 3 different studies
     c("reportedTrait", "pubmedId")]

#Investigate further...
###those that come from the same study (1:5)
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[1],c("reportedTrait")]),] #different sample size- select that one bigger N
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[2],c("reportedTrait")]),] #different sample size- select that one bigger N
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[3],c("reportedTrait")]),] #different sample size- select that one bigger N
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[4],c("reportedTrait")]),] #different sample size- select that one bigger N
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[5],c("reportedTrait")]),] #different sample size- select that one bigger N

###those that come from different studies (6:10), so we keep them at this stage
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[6],c("reportedTrait")]),] #different studies, different sample size
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[7],c("reportedTrait")]),] #different studies, different sample size and ancestry
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[8],c("reportedTrait")]),] #different studies, different sample size and ancestry
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[9],c("reportedTrait")]),] #different studies, different sample size and ancestry
GWAS[which(GWAS$reportedTrait==GWAS[which(duplicated(GWAS$reportedTrait))[10],c("reportedTrait")]),] #different studies, different sample size



#########################
#Selection of EAS GWAS
#########################

#Keep those GWAS studying participants of "East Asian" ancestry
GWAS_microbiome_EastAsians <- subset(GWAS, grepl("East Asian", discoverySampleAncestry))
table(GWAS_microbiome_EastAsians$discoverySampleAncestry, useNA = "always")

table(GWAS_microbiome_EastAsians$pubmedId, useNA = "always")

#Check on the duplicated GWAS
length(which(duplicated(GWAS_microbiome_EastAsians$reportedTrait)))
# [1] 6 gut microbiome measures appear more than once
GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait)),c("reportedTrait", "pubmedId")]
#                                                                                      reportedTrait pubmedId
# 441                                                   Gut microbiota alpha diversity (Chao1 index) 33208821
# 466                                                 Gut microbiota relative abundance (Prevotella) 33208821
# 467                                           Gut microbiota relative abundance (Faecalibacterium) 33208821
# 468                                               Gut microbiota relative abundance (Oscillospira) 33208821
# 469 Gut microbiota relative abundance (unclassified genus belonging to family Erysipelotrichaceae) 33208821
# 677                                                 Gut microbiota alpha diversity (Shannon index) 33462485
#2 studies (1 of the 6 duplicated traits must be kept at this stage because duplicated traits come from different studies, the other 5 are actually repeated for some reason - which below)
                        
table(GWAS_microbiome_EastAsians$pubmedId, GWAS_microbiome_EastAsians$discoverySampleAncestry, useNA = "always")
table(GWAS_microbiome_EastAsians$pubmedId, GWAS_microbiome_EastAsians$initialSampleDescription, useNA = "always")

GWAS_microbiome_EastAsians$discoverySampleAncestry <- ifelse(GWAS_microbiome_EastAsians$discoverySampleAncestry=="14306 European,811 East Asian,1097 Hispanic or Latin American,114 African American or Afro-Caribbean,481 Greater Middle Eastern (Middle Eastern, North African or Persian),1531 NR, Other admixed ancestry",
                                                             "811 East Asian (and other pop)",
                                                             GWAS_microbiome_EastAsians$discoverySampleAncestry)
table(GWAS_microbiome_EastAsians$discoverySampleAncestry, useNA = "always")

#Create a variable for sample size (deleting the non-numeric characters)
GWAS_microbiome_EastAsians$sample_size <- as.numeric(gsub("\\D", "", GWAS_microbiome_EastAsians$discoverySampleAncestry))
table(GWAS_microbiome_EastAsians$sample_size)

GWAS_microbiome_EastAsians$sample_size <- ifelse(
  !is.na(GWAS_microbiome_EastAsians$replicationSampleAncestry), 
  as.numeric(sub(".*?(\\d+) East Asian.*", "\\1", GWAS_microbiome_EastAsians$discoverySampleAncestry)) + as.numeric(sub(".*?(\\d+) East Asian.*", "\\1", GWAS_microbiome_EastAsians$replicationSampleAncestry)), 
  as.numeric(sub(".*?(\\d+) East Asian.*", "\\1", GWAS_microbiome_EastAsians$discoverySampleAncestry))
)
#Order by sample size
head(GWAS_microbiome_EastAsians)
tail(GWAS_microbiome_EastAsians)
GWAS_microbiome_EastAsians <- GWAS_microbiome_EastAsians[order(GWAS_microbiome_EastAsians$sample_size, decreasing = T),]
head(GWAS_microbiome_EastAsians)
tail(GWAS_microbiome_EastAsians)

#Do something about duplicated studies 
length(which(duplicated(GWAS_microbiome_EastAsians$reportedTrait)))
GWAS_microbiome_EastAsians[which(GWAS_microbiome_EastAsians$reportedTrait==GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[1],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]  #this trait comes from 2 different studies, so keep both at this stage
GWAS_microbiome_EastAsians[which(GWAS_microbiome_EastAsians$reportedTrait==GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[2],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]  #this trait comes from the same study- select that one bigger N
GWAS_microbiome_EastAsians[which(GWAS_microbiome_EastAsians$reportedTrait==GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[3],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]  #this trait comes from the same study- select that one bigger N
GWAS_microbiome_EastAsians[which(GWAS_microbiome_EastAsians$reportedTrait==GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[4],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]  #this trait comes from the same study- select that one bigger N
GWAS_microbiome_EastAsians[which(GWAS_microbiome_EastAsians$reportedTrait==GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[5],c("reportedTrait")]), #this trait comes from 2 different studies
     c("reportedTrait", "pubmedId")]  #this trait comes from the same study- select that one bigger N
GWAS_microbiome_EastAsians[which(GWAS_microbiome_EastAsians$reportedTrait==GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[6],c("reportedTrait")]), #this trait comes from 3 different studies
     c("reportedTrait", "pubmedId")] #this trait comes from the same study- select that one bigger N 

#Just keep one of the duplicated pairs (that with bigger N - dataset rows ordered by N)
###First create one subset with those unique 
GWAS_EAS_nodup <- GWAS_microbiome_EastAsians[!duplicated(GWAS_microbiome_EastAsians$reportedTrait),] 
###Then create one subset with that duplicated but from different studies (we want to keep both at this stage)
trait_dup_2gwas <- GWAS_microbiome_EastAsians[which(GWAS_microbiome_EastAsians$reportedTrait==
                                                      GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[1],c("reportedTrait")]),
                                              ]  #this trait comes from 2 different studies, so keep both at this stage
###Create a single dataset with those two above
GWAS_EAS_nodup <- rbind(GWAS_EAS_nodup,trait_dup_2gwas)
###Then create one subset with those duplicated to select only the one with bigger N
GWAS_EAS_dup <- GWAS_microbiome_EastAsians[which(GWAS_microbiome_EastAsians$reportedTrait==
                                                   GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[2],c("reportedTrait")]|
                                                   GWAS_microbiome_EastAsians$reportedTrait==
                                                   GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[3],c("reportedTrait")]|
                                                   GWAS_microbiome_EastAsians$reportedTrait==
                                                   GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[4],c("reportedTrait")]|
                                                   GWAS_microbiome_EastAsians$reportedTrait==
                                                   GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[5],c("reportedTrait")]|
                                                   GWAS_microbiome_EastAsians$reportedTrait==
                                                   GWAS_microbiome_EastAsians[which(duplicated(GWAS_microbiome_EastAsians$reportedTrait))[6],c("reportedTrait")]),
                                           ]  #this trait comes from 2 different studies, so keep both at this stage
###Check they are ordered by sample size and select the one of each pair with bigger N
GWAS_EAS_dup[,"sample_size"]
GWAS_EAS_dup <- GWAS_EAS_dup[-which(duplicated(GWAS_EAS_dup$reportedTrait)),]
GWAS_EAS_dup[,"sample_size"]
###Create a single dataset with those originally not duplicated and those selected
GWAS_EAS <- rbind(GWAS_EAS_nodup,GWAS_EAS_dup)
GWAS_EAS <- GWAS_EAS[order(GWAS_EAS$sample_size, decreasing = T),]

#Save dataset
GWAS_EAS <- GWAS_EAS[,c("firstAuthor","publicationDate","journal","title","pubmedId",                  
                        "reportedTrait",
                        "discoverySampleAncestry","initialSampleDescription","sample_size",
                        "replicationSampleAncestry","replicateSampleDescription")]
write.table(GWAS_EAS,"East_Asians_microbiome_GWAS_noduplis.csv", sep = "\t", row.names = F) 

#Save IDs of the selected GWAS
GWAS_EAS_ids <- GWAS_EAS[!duplicated(GWAS_EAS$pubmedId),c("firstAuthor","publicationDate","journal","title","pubmedId")]
write.table(GWAS_EAS_ids,"East_Asians_microbiome_GWAS_noduplis_studies.csv", sep = "\t", row.names = F) 

#Save exposures of the selected GWAS 
length(levels(as.factor(GWAS_EAS$reportedTrait)))
length(unique(levels(as.factor(GWAS_EAS$reportedTrait))))
EAS_traits <- GWAS_EAS[,c("reportedTrait")]
write.table(EAS_traits,"East_Asians_microbiome_traits.csv", sep = "\t", row.names = F, col.names = F) 
