rm(list=ls())

library(readr)
library(ggforestplot)
library(tidyverse)
library(ggplot2)
library(ggforce)

main_res <- read.csv2("../../results/MR/main_analysis_Nalls.csv") #Run this is the object is not already created

main_res_sign <- read.csv2("../../results/MR/main_analysis_Nalls_significative.csv", check.names = FALSE) #Run this as information has been added
main_res_sign <- main_res_sign[which(main_res_sign$or < 5),]
main_res_sign <- main_res_sign[which(main_res_sign$or > 0.1),]
main_res_sign <- main_res_sign %>%
    arrange(`Common class`)

#Significative results
ggforestplot::forestplot(
  df = main_res_sign,
  name = exposure,
  estimate = b,
  se = se,
  pvalue = pval,
  colour = `Common class`,
  shape = Method,
  psignif = 0.05,
  logodds = TRUE
) +
  ggforce::facet_col(
    facets = ~Gram, #Interchangable between Gram and Class or whatever characteristic
    scales = "free_y",
    space = "free"
  ) + 
  ggtitle("Bacterial significative associations to Parkinson's Disease risk") +
  xlab("Odds ratio (95% CI)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank()
  ) +
  ggplot2::scale_shape_ordinal(
    labels = c("Inverse Variance Weighted", "Wald Ratio")
  )



###MEDIATION###
il17_res <- read.csv2("mediation_IL17_PD_Nalls.csv")
b12_res <- read.csv2("mediation_B12_PD_Nalls.csv")

mediation_res <- rbind(b12_res,il17_res)

write.csv2(mediation_res,"mediation_Nalls.csv", row.names = FALSE)

mediation_res <- read.csv2("mediation_Nalls.csv", check.names = FALSE) #This file has added mediator features for plotting
mediation_res <- mediation_res %>%
  arrange(se)
ggforestplot::forestplot(
  df = mediation_res,
  name = exposure,
  estimate = b,
  se = se,
  pvalue = pval,
  colour = Biomolecule,,
  shape = Method,
#  psignif = 0.05,
  logodds = TRUE
) +
  ggforce::facet_col(
    facets = ~Mediator, #Interchangable between Gram and Class or whatever characteristic
    scales = "free_y",
    space = "free"
  ) + 
  ggtitle("Mediators associations to Parkinson's Disease risk") +
  xlab("Odds ratio (95% CI)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    legend.position = "bottom"
  ) +
  ggplot2::scale_shape_ordinal(
    labels = c("Inverse Variance Weighted", "Wald Ratio")
  )

