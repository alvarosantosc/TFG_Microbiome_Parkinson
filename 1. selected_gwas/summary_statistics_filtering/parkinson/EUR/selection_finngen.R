rm(list = ls())


library(readr)
library(dplyr)
library(vroom)


parkinson <- data.frame(vroom("summary_stats_finngen_R10_G6_PARKINSON.gz")) #Downloaded summary statistics from Finngen Parkinson Disease
head(parkinson)
names(parkinson)
str(parkinson)

#Obtaining variants in the genome wide significance threshold
parkinson_fil <- parkinson[which(parkinson$pval < 5e-8),]

#Checking duplicates
cat(c("There are",length(which(duplicated(parkinson_fil$rsids))),"duplicated rsIDs"))

parkinson_fil[which(duplicated(parkinson_fil$rsids)),]
parkinson_fil[which(parkinson_fil$rsids == "rs2696531"),]
parkinson_fil[which(parkinson_fil$rsids == "rs2696530"),]

#These are not really duplicates but SNPs in the same position whose effect allele differs, but as the clumping step does not consider this, I will keep the most relevant one

#Filtering by p-val only if this differs by more than 1 magnitude order
parkinson_fil_p <- as.data.frame(parkinson_fil %>%
                                group_by(rsids) %>%
                                mutate(highest_p_value = min(pval)) %>%
                                filter(pval / highest_p_value < 10) %>%
                                select(-highest_p_value) %>%
                                ungroup())

cat(c(dim(parkinson_fil)[1] - dim(parkinson_fil_p)[1],"duplicates have been deleted"))
#Filtering by beta value
parkinson_fil_b <- as.data.frame(parkinson_fil_p %>%
                                   group_by(rsids) %>%
                                   arrange(desc(abs(beta))) %>%
                                   slice(1) %>%
                                   ungroup())

#Checking if it worked correctly
parkinson_fil_b[which(parkinson_fil_b$rsids == "rs2696531"),"beta"]
parkinson_fil[which(parkinson_fil$rsids == "rs2696531"),"beta"]

parkinson_fil_b[which(parkinson_fil_b$rsids == "rs2696530"),"beta"]
parkinson_fil[which(parkinson_fil$rsids == "rs2696530"),"beta"]

#It did keep the highest absolute beta value

#I have identified several rsIDs that are not SNPs but indels, keeping only SNPs
parkinson_final <- parkinson_fil_b %>%
  filter(nchar(alt) == 1) %>%
  filter(nchar(ref) == 1)

write.csv(parkinson_final,"../../3.Selection_SNPs_outcome/2a.GWASs_SNPs_EU/finngen_outcome.csv", row.names = FALSE)
