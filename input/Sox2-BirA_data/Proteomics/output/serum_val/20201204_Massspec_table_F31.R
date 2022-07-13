# 20201202
# ASM 
#F31 

# BirA/Sox2Cre: Figure tables for mass spec data, F31

# last updated: 2021-02-18

################################################ read in data ##########################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(gt)

ASM_serum <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/TR01_Serum_AmandaPlex_Multimedian_notNormalized_Unfiltered_correctOrder.csv")
ASM_serum1 <- ASM_serum[grep("C7|Fstl1|Apoa1|Cfd|Adipoq|Pzp|A2m|Alb", ASM_serum$geneSymbol), ]
head(ASM_serum1[1:5, 1:5])
ASM_serum1$geneSymbol

ASM_serumScore <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/Amanda_not_normalized_filtered_20201201.csv")
ASM_serum1Score <- ASM_serumScore[grep("C7|Fstl1|Apoa1|Cfd|Adipoq|Pzp|A2m|Alb", ASM_serumScore$Gene_Symbol), ]
head(ASM_serum1Score[1:5, 1:5])
ASM_serum1Score$Gene_Symbol

write.csv(ASM_serum1, file = "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/20210218_LCMSMS_valHits.csv")
write.csv(ASM_serum1Score, file = "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/20210218_LCMSMS_valHits_MSScore_Claire.csv")

valHits <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/20210218_LCMSMS_valHits.csv", 
                    nrows = 9)
msScore_hits <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/20210218_LCMSMS_valHits_MSScore_Claire.csv", 
                         nrow =9)


msScore_hits <- msScore_hits[ , c("Gene_Symbol", "PC_NC", "Score..FPR10_20201130.")]
colnames(msScore_hits) <- c("Gene", "PC_NC", "FPR_10")
#msScore_hits$Gene <- sort(msScore_hits$Gene)
msScore_hits <- msScore_hits %>% arrange(Gene)

# fix Cfd, is listed as a PC
msScore_hits$PC_NC[3] <- "PC"

# remove accidental Albfm1 protein
msScore_hits <- msScore_hits[c(1:3, 5:9), ]

valHits$logFC.xBirA.Cre.over.NoBirA
valHits$adj.P.Val.xBirA.Cre.over.NoBirA

valHits <- valHits[ , c("id","logFC.xBirA.Cre.over.NoBirA","adj.P.Val.xBirA.Cre.over.NoBirA", "geneSymbol", "accession_number")]
valHits <- data.frame(Accession = valHits$id, Gene = valHits$geneSymbol, logFC = valHits$logFC.xBirA.Cre.over.NoBirA, adjP = valHits$adj.P.Val.xBirA.Cre.over.NoBirA)
str(valHits)

#valHits$logFC <- abs(valHits$logFC) # not needed with correctOrder data from Kiki
valHits
valHits <- valHits %>% arrange(Gene)
valHits <- cbind(valHits, msScore_hits[ , 2:3])

# re-name Cfd as Adipsin
valHits <- valHits %>% arrange(Gene)


# make gt table 
valHits_gt <- gt(valHits[ , 2:6])


# format Gene col to be bold 
valHits_gt <- valHits_gt %>%
  tab_style(
    style = cell_text(size = px(15), font = "arial"),
    locations = cells_body(vars(`geneSymbol`))
  )

valHits_gt

# change all cols format, except Gene
valHits_gt <- valHits_gt %>%
  tab_style(
    style = cell_text(
      size = px(15),
      color = "black",
      font = "arial"
    ),
    locations = cells_body(vars(`logFC`, `adjP`, `PC_NC`, `FPR_10`))
  ) 

valHits_gt

# format col headers to be all caps 
valHits_gt <- valHits_gt %>%
  tab_style(
    style = cell_text(
      size = px(15),
      color = "black",
      font = "arial"
    ),
    locations = cells_column_labels(everything())) %>% 
  tab_style(
    style = cell_text(
      size = px(15),
      color = "black",
      font = "arial"
    ),
    locations = cells_body(vars(Gene))
  ) %>% tab_style(
    style = cell_text(
      size = px(15),
      color = "black",
      font = "arial"
    ),
    locations = cells_column_labels(vars(Gene))
  ) 



valHits_gt

valHits_gt <- valHits_gt %>%
  tab_options(
    column_labels.border.top.style = "#334422",
    column_labels.border.top.width = 1,
    table.border.top.style = "none",
    column_labels.border.bottom.style = "none",
    column_labels.border.bottom.width = 1,
    column_labels.border.bottom.color = "#334422",
    table_body.border.top.style = "none",
    table_body.border.bottom.color = "#0000001A",
    data_row.padding = px(7)
  )
valHits_gt

# modify table cols to have specific widths 
valHits_gt <- valHits_gt %>%
  cols_width(
    vars(Gene) ~ px(75),
    vars(`logFC`) ~ px(110),
    vars(`adjP`) ~ px(125), 
    vars(`PC_NC`) ~ px(75), 
    vars(`FPR_10`) ~ px(100)
  )

valHits_gt

# re-label col labels, center col labels, and add tab spanners for tissue type
valHits_gt <- valHits_gt %>% 
  cols_label(logFC = html("logFC <br>BirA-Sox2Cre<br>/BirA"), adjP = "Adj. p-value", 
             PC_NC = html("Secreted<br>Postive<br>Control"), FPR_10 = html("Enrichment<br>Score"))  %>% 
  cols_align(align = "center", columns = vars(Gene, `logFC`, `adjP`, `PC_NC`, `FPR_10`)) %>% 
  tab_style(
    style = cell_text(size = px(15), weight = "bold", font = "arial"),
    locations = cells_column_labels(vars(Gene, logFC, adjP, PC_NC, FPR_10))) %>% 
  tab_footnote(
    footnote = "PC = positive control; Uniprot annotated secreted protein",
    locations = cells_column_labels(
      columns = vars(`PC_NC`))
)


valHits_gt

library(webshot)
webshot::install_phantomjs()
gtsave(valHits_gt, "20210218_BirASox2Cre_LCMSMS_WB_validatedHits_table_logFC_adjPvalue1__MSScorev2.pdf", 
       path = "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/")









