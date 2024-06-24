rm(list = ls())


library(readr)
library(dplyr)
library(readxl)
library(openxlsx)

GWAS_1 <- read.csv2("GWAS 1 microbioma EU_raw.csv", skip = 2)
head(GWAS_1)
names(GWAS_1)

#Obtaining only those with a p < 5e-8 and with information about its position
GWAS_1_p <- GWAS_1[which(GWAS_1$META.P <= 5e-8),]
GWAS_1_p <- GWAS_1_p[which(is.na(GWAS_1_p$pos) == FALSE),]

#Deleting not relevant columns
GWAS_1_p <- GWAS_1_p[,1:12]
names(GWAS_1_p)

#Checking for duplicates
length(which(duplicated(GWAS_1_p$pos)))
#14 duplicates
dups_1 <- GWAS_1_p$pos[which(duplicated(GWAS_1_p$pos))]
GWAS_1_dup <- GWAS_1_p[which(GWAS_1_p$pos %in% dups_1),]

#Getting only the duplicate with the highest absolute beta value among duplicates
GWAS_1_fil <- as.data.frame(GWAS_1_p %>%
  group_by(pos) %>%
  arrange(desc(abs(META.BETA))) %>%
  slice(1) %>%
  ungroup())

#Checking a duplicated one randomly to see if the higher beta values have been selected
GWAS_1_p[which(GWAS_1_p$pos == "239148593"),c("pos","META.BETA")]
GWAS_1_fil[which(GWAS_1_fil$pos == "239148593"),c("pos","META.BETA")]
#Yes, the code worked as intented and the highest beta values among duplicates have been selected

write.csv2(GWAS_1_fil,"GWAS 1 microbioma filtrado.csv")




GWAS_2 <- read.csv2("GWAS 2 microbioma EU_raw.csv", skip = 2)
head(GWAS_2)
names(GWAS_2)

#Although supposedly all of these SNPs are genome-wide significant, we check it
GWAS_2_p <- GWAS_2[which(GWAS_2$P <= 5e-8),]
#Actually 1 disappeared

GWAS_2[which(GWAS_2$P >= 5e-8),c("Associated.variant..rsID.","P")]
#It is weird as the df is ordered from lowest to highest P value, and this one just appears among the e-08, will not consider it

length(which(duplicated(GWAS_2_p$Associated.variant..rsID.)))
#116 duplicated SNPs
dups_2 <- GWAS_2_p$Associated.variant..rsID.[which(duplicated(GWAS_2_p$Associated.variant..rsID.))]
GWAS_2_dup <- GWAS_2_p[which(GWAS_2_p$Associated.variant..rsID. %in% dups_2),]

#Here I group all rows with the same SNP and create a new variable which is the lowest P value among each group. If that row's P value differ from the best one by more than 1 magnitude unit, it is deleted
GWAS_2_fil_p <- as.data.frame(GWAS_2_p %>%
  group_by(Associated.variant..rsID.) %>%
  arrange(Associated.variant..rsID.) %>%
  mutate(highest_p_value = min(P)) %>%
  filter(P / highest_p_value < 10) %>%
  ungroup())

cat(c(dim(GWAS_2_p)[1] - dim(GWAS_2_fil_p)[1],"duplicated SNPs have been deleted"))
#Checking if it worked
GWAS_2_fil_p[which(GWAS_2_fil_p$P / GWAS_2_fil_p$highest_p_value >= 10),"Associated.variant..rsID."]
#It did! Now we have to filter the rest of the SNPs by beta
#Deleting the added variable
GWAS_2_fil_p <- GWAS_2_fil_p[,-ncol(GWAS_2_fil_p)]

#First group by SNP, then arrange each group from highest absolute beta value to lowest, finally keep only the first one from each group
GWAS_2_fil_b <- as.data.frame(GWAS_2_fil_p %>%
  group_by(Associated.variant..rsID.) %>%
  arrange(desc(abs(beta))) %>%
  slice(1) %>%
  ungroup())

length(which(duplicated(GWAS_2_fil_b$Associated.variant..rsID.)))
#No duplicates remain

write.csv2(GWAS_2_fil_b,"GWAS 2 microbioma filtrado.csv")




GWAS_3 <- read.csv2("GWAS 3 microbioma EU_raw.csv",skip = 4)
head(GWAS_3)
names(GWAS_3)

GWAS_3 <- GWAS_3[,1:20]

GWAS_3_p <- GWAS_3[which(GWAS_3$meta_P <= 5e-8),]

#HB ones refer to presence/absence while RNT refers to abundance, so we will only keep those
GWAS_3_p <- GWAS_3_p[which(GWAS_3_p$model == "RNT"),]

length(which(duplicated(GWAS_3_p$rsid)))
#When filtering for the genome-wide associated SNPs no duplicates remain, so that's all

write.csv2(GWAS_3_p,"GWAS 3 microbioma filtrado.csv")




GWAS_4 <- as.data.frame(read_excel("GWAS 4 microbioma EU_raw.xlsx"))
head(GWAS_4)
names(GWAS_4)

convert_to_scientific <- function(x) {
  # Replace special characters with standard notation
  x <- gsub(" × ", "e", x)  # Replace " × " with "e"
  x <- gsub("10−", "-", x) # Replace "−" with "-"
  return(x)
}

GWAS_4$`Meta P` <- sapply(GWAS_4$`Meta P`, convert_to_scientific)

GWAS_4$`Meta P` <- as.numeric(GWAS_4$`Meta P`)

GWAS_4_p <- GWAS_4[GWAS_4$`Meta P` <= 5e-8,]

length(which(duplicated(GWAS_4_p$SNP)))
#12 duplicated SNPs

GWAS_4_p[which(duplicated(GWAS_4_p$SNP)),"SNP"]

#There is an issue with this GWAS that duplicated SNPs just show information about location in one of the duplicaates

GWAS_4_p <- GWAS_4_p %>%
  group_by(SNP) %>%
  mutate(
    Locus = first(Locus),
    Chr. = first(Chr.),
    `Locus start` = first(`Locus start`),
    `Locus end` = first(`Locus end`),
    `Nearest gene` = first(`Nearest gene`),
    `Genes in locus` = first(`Genes in locus`)) %>%
  ungroup()

#Filtering for the duplicates with a P value higher than 1 magnitude order
GWAS_4_fil_p <- as.data.frame(GWAS_4_p %>%
                                group_by(SNP) %>%
                                mutate(highest_p_value = min(`Meta P`)) %>%
                                filter(`Meta P`/ highest_p_value < 10) %>%
                                ungroup())

cat(c(dim(GWAS_4_p)[1] - dim(GWAS_4_fil_p)[1],
      "duplicated SNPs have been deleted,",
      length(which(duplicated(GWAS_4_fil_p$SNP))),
      "duplicates remain"))

#Apparently the whole dataset used lines as the minus symbol
GWAS_4_fil_p$`Meta β` <- gsub("−","-",GWAS_4_fil_p$`Meta β`)
GWAS_4_fil_p$`Meta β` <- as.numeric(GWAS_4_fil_p$`Meta β`)
GWAS_4_fil_b <- as.data.frame(GWAS_4_fil_p %>%
                                group_by(SNP) %>%
                                arrange(desc(abs(`Meta β`))) %>%
                                slice(1) %>%
                                ungroup())

cat(c(dim(GWAS_4_fil_p)[1] - dim(GWAS_4_fil_b)[1],
      "duplicated SNPs have been deleted,",
      length(which(duplicated(GWAS_4_fil_b$SNP))),
      "duplicates remain"))

GWAS_4_fil_b <- GWAS_4_fil_b[,-ncol(GWAS_4_fil_p)]

#I am not sure if the P/beta column is relevant, but as it was not correctly calculated due to the uncorrect format of the columns I will recalculate it
GWAS_4_fil_b$`β-div P` <- GWAS_4_fil_b$`Meta β` / GWAS_4_fil_b$`Meta P`

#Renaming certain columns that in csv are problematic
names(GWAS_4_fil_b)[6] <- "Meta P"
names(GWAS_4_fil_b)[7] <- "Meta Beta"
names(GWAS_4_fil_b)[8] <- "Beta div P"

write.csv2(GWAS_4_fil_b,"GWAS 4 microbioma EU.csv")


GWAS_pool <- read.xlsx("Pool de GWAS microbioma EU.xlsx")
head(GWAS_pool)
names(GWAS_pool)


length(which((duplicated(GWAS_pool$rsID) & is.na(GWAS_pool$rsID) == FALSE) |
              ((duplicated(GWAS_pool$chromosome) & is.na(GWAS_pool$chromosome) == FALSE) &
                 duplicated(GWAS_pool$position) & is.na(GWAS_pool$position) == FALSE))
)
