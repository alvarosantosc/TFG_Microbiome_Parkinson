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

GWAS <- read_tsv("GWAS microbiome_raw.tsv") #Downloaded studies from GWAS Catalog for EFO_0007874
GWAS <- as.data.frame(GWAS)
head(GWAS)
names(GWAS)
dim(GWAS)
dim(GWAS)


#########################
#QC 
#########################

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
#Selection of EUR GWAS
#########################

#Keep those GWAS studying participants of "European" ancestry
GWAS_EU <- subset(GWAS, grepl("European", discoverySampleAncestry))
table(GWAS_EU$discoverySampleAncestry, useNA = "always")
table(GWAS_EU$pubmedId, useNA = "always")

#Check on the duplicated GWAS
length(which(duplicated(GWAS_EU$reportedTrait)))
# [1] 3 gut microbiome measures appear more than once
GWAS_EU[which(duplicated(GWAS_EU$reportedTrait)),c("reportedTrait", "pubmedId")]
# reportedTrait pubmedId
# 471  Gut microbiota (bacterial taxa) 27694959
# 472  Gut microbiota (beta diversity) 27723756
# 2304 Gut microbiota (bacterial taxa) 27723756
#2 studies (1 of the 3 duplicated traits must be kept at this stage because duplicated traits come from different studies, the other 2 are actually repeated for some reason - which below)

GWAS_EU$discoverySampleAncestry <- ifelse(GWAS_EU$discoverySampleAncestry=="14306 European,811 East Asian,1097 Hispanic or Latin American,114 African American or Afro-Caribbean,481 Greater Middle Eastern (Middle Eastern, North African or Persian),1531 NR, Other admixed ancestry",
                                                             "14306 European (and other pop)",
                                          GWAS_EU$discoverySampleAncestry)

table(GWAS_EU$discoverySampleAncestry, useNA = "always")

#Create a variable for sample size (deleting the non-numeric characters)
GWAS_EU$sample_size <- as.numeric(gsub("\\D", "", GWAS_EU$discoverySampleAncestry))
table(GWAS_EU$sample_size)

GWAS_EU$sample_size <- ifelse(
  !is.na(GWAS_EU$replicationSampleAncestry), 
  as.numeric(gsub("^(\\d+) European.*", "\\1", GWAS_EU$discoverySampleAncestry)) + as.numeric(gsub("^(\\d+) European.*", "\\1", GWAS_EU$replicationSampleAncestry)), 
  as.numeric(gsub("^(\\d+) European.*", "\\1", GWAS_EU$discoverySampleAncestry))
)

#Order by sample size
head(GWAS_EU)
tail(GWAS_EU)
GWAS_EU <- GWAS_EU[order(GWAS_EU$sample_size, decreasing = T),]
head(GWAS_EU)
tail(GWAS_EU)

table(GWAS_EU$pubmedId, GWAS_EU$discoverySampleAncestry, useNA = "always")
table(GWAS_EU$pubmedId, GWAS_EU$initialSampleDescription, useNA = "always")

#What to do about duplicated studies 
length(which(duplicated(GWAS_EU$reportedTrait)))
GWAS_EU[which(GWAS_EU$reportedTrait==GWAS_EU[which(duplicated(GWAS_EU$reportedTrait))[1],c("reportedTrait")]), #this trait comes from 2 different studies
                           c("reportedTrait", "pubmedId")]  #this trait comes from 3 different studies, so keep all at this stage
GWAS_EU[which(GWAS_EU$reportedTrait==GWAS_EU[which(duplicated(GWAS_EU$reportedTrait))[2],c("reportedTrait")]), #this trait comes from 2 different studies
                           c("reportedTrait", "pubmedId")] #this trait comes from 2 different studies, so keep all at this stage
GWAS_EU[which(GWAS_EU$reportedTrait==GWAS_EU[which(duplicated(GWAS_EU$reportedTrait))[3],c("reportedTrait")]), #this trait comes from 2 different studies
                           c("reportedTrait", "pubmedId")]  #this trait comes from 3 different studies, so keep all at this stage

length(unique(GWAS_EU$pubmedId))
table(GWAS_EU$pubmedId, useNA = "always")

#Save dataset
GWAS_EU <- GWAS_EU[,c("firstAuthor","publicationDate","journal","title","pubmedId",                  
                        "reportedTrait",
                        "discoverySampleAncestry","initialSampleDescription","sample_size",
                        "replicationSampleAncestry","replicateSampleDescription")]
write.table(GWAS_EU,"Europeans_microbiome_GWAS_noduplis.csv", sep = "\t", row.names = F) 

#Get only those studies with a sample size higher than 1000
GWAS_EU_minN1000 <- as.data.frame(subset(GWAS_EU, sample_size >= 1000))
table(GWAS_EU_minN1000$sample_size)
length(unique(GWAS_EU_minN1000$pubmedId))
table(GWAS_EU_minN1000$pubmedId, useNA = "always")

write.table(GWAS_EU_minN1000,"Europeans_microbiome_GWAS_noduplis_N1000.csv", sep = "\t", row.names = F) 

#Save PubmedIDs of the selected GWAS
GWAS_EU_minN1000_ids <- GWAS_EU_minN1000[!duplicated(GWAS_EU_minN1000$pubmedId),c("firstAuthor","publicationDate","journal","title","pubmedId")]
write.table(GWAS_EU_minN1000_ids,"Europeans_microbiome_GWAS_noduplis_studies.csv", sep = "\t", row.names = F) 

#Save exposures of the selected GWAS 
length(levels(as.factor(GWAS_EU_minN1000$reportedTrait)))
length(unique(levels(as.factor(GWAS_EU_minN1000$reportedTrait))))
EUR_traits <- GWAS_EU_minN1000[,c("reportedTrait")]
write.table(EUR_traits,"Europeans_microbiome_traits.csv", sep = "\t", row.names = F, col.names = F) 
