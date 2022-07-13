# 2021-10-15
# TR01: BirA/Sox2 Serum 
# ASM 

# usage: script to calculate MS-Score (from Perrimon Group (Claire Hu)) for a 6plex LC-MSMS run of quantitative tandem mass spectrometry
# from the allRatios files from the Carr lab 



####################################################################### begin ##########################################################################
library(dplyr)
library(tidyr)
library(tidyverse)

# data from two separate runs: sample prep by Ilia (I) or Amanda (A)


# will start with ASM data to make pipeline, then will run same analysis for Ilia's data 
#serum_A <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/TR01_Serum_AmandaPlex_AllRatios_notNormalized_Unfiltered.csv") # serum_A 
serum_A <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/TR01_Serum_AmandaPlex_AllRatios_notNormalized_Score_greater25.csv") # serum_A 
head(serum_A)

# select out BirA+/BirA- cols (total 9 ratios from each BirA+ (n = 3) over each BirA- (n = 3))
colnames(serum_A)
ratio_cols <- c("Amanda_newDatabase..131.128..BirA.Cre_1.NoBirA_4", "Amanda_newDatabase..131.129..BirA.Cre_1.NoBirA_5", "Amanda_newDatabase..131.126..BirA.Cre_1.NoBirA_6", 
                "Amanda_newDatabase..130.128..BirA.Cre_2.NoBirA_4", "Amanda_newDatabase..130.129..BirA.Cre_2.NoBirA_5", "Amanda_newDatabase..130.126..BirA.Cre_2.NoBirA_6", 
                "Amanda_newDatabase..127.128..BirA.Cre_3.NoBirA_4", "Amanda_newDatabase..127.129..BirA.Cre_3.NoBirA_5", "Amanda_newDatabase..127.126..BirA.Cre_3.NoBirA_6", 
                "id", "accession_number", "geneSymbol", "entry_name", "cellular_component")
serum_A <- serum_A[ , ratio_cols]
head(serum_A[ , 1:9])


# serum BirA+/BirA- log2FC data
#serum_logFC <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/TR01_Serum_AmandaPlex_Multimedian_notNormalized_Unfiltered_correctOrder.csv")
serum_logFC_all <- read.csv("input/TR01_Plasma_AllDataTogether_Unfiltered_Unnormalized_Multimedian_CorrectDatabase.csv")
# impute NAs to 0.000
serum_A[is.na(serum_A)] <- 0
serum_logFC[is.na(serum_logFC)] <- 0.001

serum_cols <- colnames(serum_A)

## read in UniProt datasets for positive secreted protein controls (PC) and negative controls (NC)
uniP_secreted <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Protein_localization_controls/uniprot-keyword__Secreted+(KW-0964)_-filtered-organism__Mus+muscul-- copy.csv")
uniP_nuclear <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Protein_localization_controls/uniprot-keyword__Nucleus+(KW-0539)_-filtered-organism__Mus+musculu-- copy.csv")
uniP_cyto <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Protein_localization_controls/uniprot-keyword__Cytoplasm+(KW-0963)_-filtered-organism__Mus+muscu-- copy.csv")
uniP_ER <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Protein_localization_controls/uniprot-keyword__Endoplasmic+reticulum+(KW-0256)_-filtered-organism_--.csv")
head(uniP_secreted)

# PC & NC lists from Claire
PC_CY <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/Mouse_PC_uniprot_secreted.csv")
NC_CY <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/Mouse_NC_combined.csv")

# PC & NC lists from Rui 
PC_RY <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/PC-NC_lists/Mouse_PC_uniprot_secreted_RY.csv")
NC_RY <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/PC-NC_lists/Mouse_NC_combined_updated_RY.csv")

############################################################################# PC/NC ID ###################################################################################
ifelse(serum_A$accession_number == PC_RY$Entry & serum_A$geneSymbol != NC_RY$Gene.Symbol, "PC", "other")
ifelse(serum_A$accession_number != PC_RY$Entry & serum_A$geneSymbol == NC_RY$Gene.Symbol, "NC", "other")

serum_A$PCorNC <- serum_A

######################################################################### Thresh-holding ##########################################################################





















serum_hits1 <- serum_hits[grep("extracellular space|extracellular exosome|extracellular", serum_hits$cellular_component), ]
str(serum_hits1)
length(serum_hits1$id)
head(serum_hits1[1:5, 1:5])
length(grep("extracellular space|extracellular exosome|extracellular", serum_hits$cellular_component))
length(grep("extracellular space|extracellular exosome|extracellular", serum_hits$cellular_component))

PCs <- intersect(PC_CY$Gene.names, serum_A$geneSymbol)
NCs <- intersect(NC_CY$Gene.Symbol, serum_A$geneSymbol)
thresh <- length(PCs)/(length(PCs) + length(NCs))

serum_131.126 <- serum_A[ , c("Amanda_newDatabase..131.126..BirA.Cre_1.NoBirA_6", "geneSymbol")]
serum_131.126$CTR <- any[]
serum_131.126$thresh <- 
  for (i in serum_131.126$Amanda_newDatabase..131.126..BirA.Cre_1.NoBirA_6) { 
  thresh_1 <- thresh/i
  }
  
thresh_df <- data.frame()
for (i in (serum_131.126$Amanda_newDatabase..131.126..BirA.Cre_1.NoBirA_6)) { 
  thresh_1 <- thresh/i
  print(i)
  print(thresh_1)
  #rbind(thresh_df, thresh_1)
  #return(thresh_df)
}

for (i in serum_131.126$Amanda_newDatabase..131.126..BirA.Cre_1.NoBirA_6) {
  
  print(i)
}

library(ggplot2)
ggplot(serum_A) + 
  geom_histogram(aes(threshold), binwidth = 0.0005, bins = 1) + 
  xlim(-0.1, 0.1)



# comparing via accession numbers
BirA_1_PC <- (length(which((serum_A$id) %in% (uniP_secreted$Entry))))/(length(uniP_secreted$Entry) + length(uniP_cyto$Entry) + length(uniP_nuclear$Entry) + length(uniP_ER$Entry))

# find number of matching PCs in TMT allRatios compared to PC list 
length(which((serum_A$id) %in% (uniP_secreted$Entry)))
(length(uniP_secreted$Entry) + length(uniP_cyto$Entry) + length(uniP_nuclear$Entry) + length(uniP_ER$Entry))

# plot data: divide % PC by each proteins logFC 
for (i in serum_logFC$logFC.xBirA.Cre.over.NoBirA) {
  plot_data <- BirA_1_PC/i
  print(plot_data)
  return(plot_data)
}

serum_logFC$logFC.xBirA.Cre.over.NoBirA
BirA_1_PC

### not working... return df incorrect lenght. Loop works outside of function 
threshold_fnct <- function(df_colname, PC_list) { 
  df10 <- data.frame()
  for (i in df_colname) {
    plot_data <- PC_list/i
    print(plot_data)
    #return(plot_data)
    df10 <- rbind(df10, plot_data)
  }
  return(df10)
}

threshold_fnct(serum_logFC$logFC.xBirA.Cre.over.NoBirA, BirA_1_PC)
df10

length(serum_logFC$logFC.xBirA.Cre.over.NoBirA)
length(df10$X0.971614974988315)

# calculate thresholds
df11 <- data.frame()
for (i in serum_logFC$logFC.xBirA.Cre.over.NoBirA) {
  plot_data <- BirA_1_PC/i
  print(i)
  print(plot_data)
  #return(plot_data)
  df11 <- rbind(df11, plot_data)
}
df11
length(df11$X0.00287522792738507)



# attach values, needs to be atttached as numeric and not a df 
serum_logFC$threshold <- as.numeric(df11$X0.00287522792738507)

str(serum_logFC)


## imouted values to 0.001 are threshold 12.206 values 
library(ggplot2)
ggplot(serum_logFC) + 
  geom_histogram(aes(threshold), binwidth = 0.0005, bins = 1) + 
  xlim(-0.1, 0.1)


ggplot(serum_logFC, aes(threshold, id)) 




?geom_histogram
?ggplot()



