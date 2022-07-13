# 2021-10-10
# TR01: BirA/Sox2 Serum, Liver, Kidney, Brain
# ASM 

# last updated: 2022-02-18

# usage: use upset plot to show proteins shared and unique to brain, liver, kidney, and serum Sox2-BirA MS datasets. 
        # used for manuscript figure 4B. 


####################################################################### begin ##########################################################################
library(dplyr)
library(tidyr)
library(tidyverse)

# load in all datasets
serum <- read.csv("serum_enriched_combined_from_RY.csv")
brain <- read.csv("BRNdata.csv") 
brain <- brain[brain$ES_score >= 5, ] # filter for enriched proteins (ES >= 5)
kid <- read.csv("KIDdata.csv")
kid <- kid[kid$ES_score >= 5, ] # filter for enriched proteins (ES >= 5)
liv <- read.csv("LIVdata.csv")
liv <- liv[liv$ES_score >= 5, ] # filter for enriched proteins (ES >= 5)

####################################################################### upset plots: tissue of origin comp ##########################################################################
# use all enriched for tissues and multimedian for serum 
library(UpSetR)

# make input list
tissue_list <- list(Serum = serum$geneSymbol, Brain = brain$geneSymbol, Kidney = kid$geneSymbol, Liver = liv$geneSymbol) 
plot1 <- upset(fromList(tissue_list), order.by = "freq")
#str(plot1)
plot1


