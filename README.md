This repository contains all the code runned in the research project "Characterising the causal effects of the Gut Microbiota on Parkinson’s Disease risk: A Mendelian Randomization study"

#################
**Analysis overview**
#################

Mendelian Randomization Study of gut bacteria on Parkinson's disease (PD) and follow-up analyses

################add graphical abstract

***Analysis components***
#################

1.- Hypothesis-free Mendelian Randomization discovery study - gut bacterial taxa analysed in publicly available GWAS

2.- Hypothesis-driven follow up studies - sensitivity and mediation analyses

***Required Data files***
#################

- Summary statistics of the exposure GWAS (gub microbiota)
- Summary statistics of the outcome GWAS (PD)
- File identifying the bacteria to study (exposure)
- Summary statistics of the GWAS of IL-17 proteins
- Summary statistics of the GWAS of vitamin B12

***Required software and R packages***
#################

- R ###version
- TwoSampleMR
- ################add all necessary packages

#################
**Directory structure**
#################

1) Selected GWAS: This folder is divided into 2 subfolders, one for the filtering of the studies extracted from GWAS Catalog and a second one for the SNPs filtering of these studies.
2) Main_analysis: This folder contains every step of the main analysis and the sensitivity analysis on EAS population. It is organised as: exposure data preparation, outcome data preparation, harmonization and the MR analysis that lead to the results exposed.
3) Reverse MR: Contains the code for the reverse MR analysis until the instrumental variable selection step, where no IVs were obtained.
4) Mediation: Contains all the code runned to obtain the results for mediatiors.
5) Data visualization: Consists of the code runned to obtain the exposed figures.

