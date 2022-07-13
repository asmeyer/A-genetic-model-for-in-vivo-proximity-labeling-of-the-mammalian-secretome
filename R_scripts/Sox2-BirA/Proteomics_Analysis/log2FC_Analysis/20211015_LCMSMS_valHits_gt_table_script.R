# 20211015
# ASM 

# BirA/Sox2Cre: Figure tables for mass spec data

# last updated: 2022-02-05

################################################ read in data ##########################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(gt)

########## use 
All_serum <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/input/TR01_USE_Plasma_AllDataTogether_Filtered_notNormalized_Multimedian_Correct.csv")
All_serum1 <- All_serum[grep("C7|Fstl1|Apoa1|Cfd|Adipoq|Pzp|Umod|Agt|Gpx3|Oxt|Alb|Umod", All_serum$geneSymbol), ]
head(All_serum1[1:5, 1:5])
All_serum1$geneSymbol

ASM_serum <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/input/TR01_Serum_AmandaPlex_Multimedian_notNormalized_Unfiltered_correctOrder.csv")
ASM_serum1 <- ASM_serum[grep("C7|Fstl1|Apoa1|Cfd|Adipoq|Pzp|Umod|Agt|Gpx3|Oxt|Alb|Umod|A2m|Igfbp2", ASM_serum$geneSymbol), ]
head(ASM_serum1[1:5, 1:5])
ASM_serum1$geneSymbol


### use ID serum for Umod & all others
ID_serum <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/input/TR01_Serum_IliaPlex_Multimedian_notNormalized_Unfiltered_correctOrder.csv")
ID_serum1 <- ID_serum[grep("C7|Apoa1|Cfd|Adipoq|Pzp|Umod|Agt|Gpx3|Alb|Umod|A2m|Igfbp2", ID_serum$geneSymbol), ]
head(ID_serum1[1:5, 1:5])
ID_serum1$geneSymbol

ASM_serumScore <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/input/ES/Amanda_not_normalized_filtered_20201201.csv")
ASM_serum1Score <- ASM_serumScore[grep("C7|Fstl1|Apoa1|Cfd|Adipoq|Pzp|Umod|Agt|Gpx3|Oxt|Alb|Umod|A2m|Igfbp2", ASM_serumScore$Gene_Symbol), ]
head(ASM_serum1Score[1:5, 1:5])
ASM_serum1Score$Gene_Symbol

# save protein lists
write.csv(ASM_serum1, file = "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/serum_val/20211015_LCMSMS_serum_valHits_v6.csv")
write.csv(ASM_serum1Score, file = "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/serum_val/20211015_LCMSMS_serum_valHits_MSScore_Claire_v6.csv")

# save protein lists
write.csv(ID_serum1, file = "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/serum_val/20211015_LCMSMS_serum_valHits_v6_ID.csv")


#################### fix for UMOD later 
ID_serum1$geneSymbol
ID_serum1Score$Gene_Symbol

# read in protein lists 
valHits <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/serum_val/20211015_LCMSMS_serum_valHits_v6_ID.csv")
msScore_hits <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/serum_val/20211015_LCMSMS_serum_valHits_MSScore_Claire_v6.csv")

msScore_hits <- msScore_hits[ , c("Gene_Symbol", "PC_NC", "Score..FPR10_20201130.")]
colnames(msScore_hits) <- c("Gene", "PC_NC", "FPR_10")
#msScore_hits$Gene <- sort(msScore_hits$Gene)
msScore_hits <- msScore_hits %>% arrange(Gene)
msScore_hits


############ remove accidental Albfm1 protein
msScore_hits <- msScore_hits[c(1:4, 6:8, 10, 11, 13), ]
msScore_hits

# fix Cfd, is listed as a PC
msScore_hits$PC_NC[3] <- "PC"
msScore_hits$PC_NC[6] <- "PC"
#msScore_hits$PC_NC[9] <- "PC"
msScore_hits

# add in NA for Umod 
UMOD <- data.frame(Gene = "Umod", PC_NC = "PC", FPR_10 = NA)
msScore_hits <- rbind(msScore_hits, UMOD)


valHits$logFC.xBirA.Cre.over.NoBirA
valHits$adj.P.Val.xBirA.Cre.over.NoBirA
valHits$geneSymbol



# select columns 
valHits <- valHits[ , c("id","logFC.xBirA.Cre.over.NoBirA","adj.P.Val.xBirA.Cre.over.NoBirA", "geneSymbol", "accession_number")]
valHits <- data.frame(Accession = valHits$id, Gene = valHits$geneSymbol, logFC = valHits$logFC.xBirACre.over.NoBirA, adjP = valHits$adj.P.Val.xBirACre.over.NoBirA)
str(valHits)

#valHits$logFC <- abs(valHits$logFC) # not needed with correctOrder data from Kiki
valHits <- valHits %>% arrange(geneSymbol)
valHits

### combine ES and log2FC data
valHits1 <- cbind(valHits, msScore_hits[ , 2:3])
colnames(valHits1) <- c("id", "logFC", "adjP", "Gene", "accession", "PC_NC", "FPR_10")
#valHits1 <- valHits1 %>% arrange(Gene)

# re-name Cfd as Adipsin
#valHits1["Gene"]["Cfd"] <- "Adipsin"
#valHits1 <- valHits1 %>% arrange(Gene)

########################################################## gt make table ##################################################################
# make gt table 
library(gt)
valHits_gt <- gt(valHits1[ , c(2:4, 6, 7)])


# format Gene col to be bold 
valHits_gt <- valHits_gt %>%
  tab_style(
    style = cell_text(size = px(15), font = "arial"),
    locations = cells_body(c(`Gene`))
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
    locations = cells_body(c(`logFC`, `adjP`, `PC_NC`, `FPR_10`))
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
    locations = cells_body(c(Gene))
  ) %>% tab_style(
    style = cell_text(
      size = px(15),
      color = "black",
      font = "arial"
    ),
    locations = cells_column_labels(c(Gene))
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
    c(Gene) ~ px(75),
    c(`logFC`) ~ px(110),
    c(`adjP`) ~ px(125), 
    c(`PC_NC`) ~ px(75), 
    c(`FPR_10`) ~ px(100)
  )

valHits_gt

# re-label col labels, center col labels, and add tab spanners for tissue type
valHits_gt <- valHits_gt %>% 
  cols_label(logFC = html("log2 FC <br>Sox2-BirA*G3<br>/Control"), adjP = "Adj. p-value", 
             PC_NC = html("Secreted<br>Postive<br>Control"), FPR_10 = html("Enrichment<br>Score"))  %>% 
  cols_align(align = "center", columns = c(Gene, `logFC`, `adjP`, `PC_NC`, `FPR_10`)) %>% 
  tab_style(
    style = cell_text(size = px(15), weight = "bold", font = "arial"),
    locations = cells_column_labels(c(Gene, logFC, adjP, PC_NC, FPR_10))) %>% 
  tab_footnote(
    footnote = "PC = positive control; Uniprot/NCBI annotated secreted protein",
    locations = cells_column_labels(
      columns = c(`PC_NC`))
  )


valHits_gt





# re-label col labels, center col labels, and add tab spanners for tissue type
valHits_gt <- valHits_gt %>% 
  cols_label(logFC = html("log2 FC <br>Sox2-BirA/<br>Control"), adjP = "Adj. p-value")  %>% 
  cols_align(align = "center", columns = c(Gene, `logFC`, `adjP`)) %>% 
  tab_style(
    style = cell_text(size = px(15), weight = "bold", font = "arial"),
    locations = cells_column_labels(c(Gene, logFC, adjP))) 


valHits_gt

library(webshot)
webshot::install_phantomjs()
gtsave(valHits_gt, "20220128_BirASox2Cre_LCMSMS_WB_validatedHits_table_logFC_adjPvalue1_v7.pdf", 
       path = "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/serum_val/")









