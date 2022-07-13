# 2021-09-11
# ASM 
# TR01: BirA/Sox2Cre paper 

# usage: get accession numbers for all significantly enriched proteins, look at protein overlap by tissue using upset plots, and get 
# data for proteins with predicted NLS (nuclear localization sequence). 

# last updated: 2021-09-11


################################################################## load data for mulimed & ES ####################################################################
library(ggvenn)
library(ggplot2)

# read in filtered multimedian data; log2FC > 1.0 and adjusted p < 0.05
serum_mlt <- read.csv("input/TR01_USE_Plasma_AllDataTogether_Filtered_notNormalized_Multimedian_Correct.csv")
serum_mlt.ft <- serum_mlt[serum_mlt$logFC.xBirACre.over.NoBirA >= 1.0 & serum_mlt$adj.P.Val.xBirACre.over.NoBirA <= 0.05, ] # could use mKate2 as cutoff

write.csv(serum_mlt.ft, "input/filtered_multmed/serum_filtered_logFC-adjP.csv")
serum_mlt.ft <- read.csv("input/filtered_multmed/serum_filtered_logFC-adjP.csv") # not normalized, still need to filter 
serum_mlt.ft$geneSymbol <- with(serum_mlt.ft, ifelse(geneSymbol == "", accession_number, geneSymbol))

# MM data for each tissue 
brain_mlt.ft <- read.csv("input/filtered_multmed/brain_filtered_logFC-adjP.csv")
brain_mlt.ft$geneSymbol <- with(brain_mlt.ft, ifelse(geneSymbol == "", accession_number, geneSymbol))

liv_mlt.ft <- read.csv("input/filtered_multmed/liv_filtered_logFC-adjP.csv")
liv_mlt.ft$geneSymbol <- with(liv_mlt.ft, ifelse(geneSymbol == "", accession_number, geneSymbol))

kid_mlt.ft <- read.csv("input/filtered_multmed/kid_filtered_logFC-adjP.csv")
kid_mlt.ft$geneSymbol <- with(kid_mlt.ft, ifelse(geneSymbol == "", accession_number, geneSymbol))

# mult med uniques
all.mm <- read.csv("output/multmed_lists/allTissue_venn_info_wSerum.csv") 
brain.mm <- all.mm$Brain
serum.mm <- all.mm$Serum
liv.mm <- all.mm$Liver
kid.mm <- all.mm$Kidney


# ES data for each tissue (unique hits)
# brain 
brain.es <- read.csv("input/ES/BRN_unique.csv")
brain.es$geneSymbol <- with(brain.es, ifelse(geneSymbol == "", accession_number, geneSymbol))


# serum
serum.es <- read.csv("input/ES/Amanda_not_normalized_filtered_20201201.csv")
serum.es$Symbol_Original <- with(serum.es, ifelse(Symbol_Original  == "", accession, Symbol_Original))

# liver
liv.es <- read.csv("input/ES/LIV_unique.csv")
liv.es$Symbol_Original <- with(liv.es, ifelse(geneSymbol  == "", accession, geneSymbol))

# kidney 
kid.es <- read.csv("input/ES/KID_unique.csv")
kid.es$Symbol_Original <- with(kid.es, ifelse(geneSymbol  == "", accession, geneSymbol))

## ALL hits with ES >= 5 (NOT uniques)
brain.es1 <- read.csv("input/ES/BRNdata.csv")
brain.es1$geneSymbol <- with(brain.es1, ifelse(geneSymbol == "", accession_number, geneSymbol))
brain.es1 <- brain.es1[brain.es1$ES_score >= 5, ]

liv.es1 <- read.csv("input/ES/LIVdata.csv")
liv.es1$geneSymbol <- with(liv.es1, ifelse(geneSymbol  == "", accession, geneSymbol))
liv.es1 <- liv.es1[liv.es1$ES_score >= 5, ]

kid.es1 <- read.csv("input/ES/KIDdata.csv")
kid.es1$geneSymbol <- with(kid.es1, ifelse(geneSymbol  == "", accession, geneSymbol))
kid.es1 <- kid.es1[kid.es1$ES_score >= 5, ]
################################################################## all accessions ####################################################################
### get all accession numbers for all enriched proteins by tissue 
# brain
brain.acc <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Brain/data/TR01_Brain_MedianNormalized_Multimedian.csv")
test2 <- brain.acc[brain.acc$geneSymbol %in% brain_mlt.ft$geneSymbol, ]
test2$id
test3 <- brain.acc[brain.acc$geneSymbol %in% brain.es1$geneSymbol, ]
test3$id
test4 <- unique(test2$id, test3$id)
write.csv(test4, "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/brain_enriched_v1.csv")

# liver
liv.acc <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Liver/data/TR01_Liver_MedianNormalized_Multimedian.csv")
test2 <- liv.acc[liv.acc$geneSymbol %in% liv_mlt.ft$geneSymbol, ]
test2$id
test3 <- liv.acc[liv.acc$geneSymbol %in% liv.es1$geneSymbol, ]
test3$id
test4 <- unique(test2$id, test3$id)
write.csv(test4, "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/liver_enriched_v1.csv")

# kidney
kid.acc <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Kidney/data/TR01_Kidney_mediannormalized_multimedian.csv")
test2 <- kid.acc[kid.acc$geneSymbol %in% kid_mlt.ft$geneSymbol, ]
test2$id
test3 <- kid.acc[kid.acc$geneSymbol %in% kid.es1$geneSymbol, ]
test3$id
test4 <- unique(test2$id, test3$id)
write.csv(test4, "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/kidney_enriched_v1.csv")

# serum 
ser.acc <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/TR01_Serum_AmandaPlex_Multimedian_notNormalized_Unfiltered_ASM.csv")
test2 <- ser.acc[ser.acc$geneSymbol %in% serum_mlt.ft$geneSymbol, ]
test2$id
test3 <- ser.acc[ser.acc$geneSymbol %in% serum.es$Symbol_Original, ]
test3$id
test4 <- unique(test2$id, test3$id)
write.csv(test4, "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/serum_enriched_v1.csv")


#################################################### Show shared and unique sig proteins by Upset Plots #########################################################################################
# upset plots 
library(UpSetR)

# read in dara from brain, liver, and kidney all proteins with ES >= 5
# brain
brain <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Brain/data/BRNdata.csv")
brain <- brain[brain$ES_score >= 5, ]

# liver
liv <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Liver/data/LIVdata.csv")
liv <- liv[liv$ES_score >= 5, ]

# kidney 
kid <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Kidney/data/KIDdata.csv")
kid <- kid[kid$ES_score >= 5, ]

# serum
serum <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/input/filtered_multmed/serum_filtered_logFC-adjP.csv")
serum <- read.csv("input/serum_enriched_combined_from_RY.csv")

# serum
# use multimed serum data from above

# test for generating upset plot (made up data)
test_list <- list(brain = c("Astn2", "CamKII", "Reln"), serum = c("Reln", "Bche", "Alb", "C7", "C3", "Slc38a10", "Egfr", "Adipoq"), liver = c("Alb", "C7", "Calr"), 
                  kidney = c("Umod", "Calu", "Slc26a4", "Slc38a10"))

# actual data for significant protein upset plot 
test_list <- list(Brain = brain$geneSymbol, Serum = serum$geneSymbol, Liver = liv$geneSymbol, 
                  Kidney = kid$geneSymbol)

## do same thing but add in RNAseq data (Brain_RNA, Brain_Protein, Liver_RNA, Liver_Protein, etc.)
library(UpSetR)
(plot1 <- upset(fromList(test_list), order.by = "freq", text.scale = c(4.5, 4.5, 3, 3, 4.5, 4.5)))

# set.scale = c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
str(plot1)
plot1



?upset


################################################################## NLS predictions ####################################################################
### get NLS predictions from all enriched proteins per tissue
# from: https://rostlab.org/services/nlsdb/




test2 <- ser.acc[ser.acc$geneSymbol %in% serum_mlt.ft$geneSymbol, ]

# read in NLS predictions (accession numbers)
brn.nls <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/NLS_predictions/brain_all_accessions_query_NLS.csv")
kid.nls <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/NLS_predictions/kidney_all_accessions_query_NLS.csv")
liv.nls <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/NLS_predictions/liver_all_accessions_query_NLS.csv")
ser.nls <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/NLS_predictions/serum_all_accessions_query_NLS.csv")

library(dplyr)
# brain
test6 <- data.frame("Query" = brain.acc$id, "Gene" = brain.acc$geneSymbol, brain.acc$Log.P.Value.xBirA.over.NoBirA, brain.acc$adj.P.Val.xBirA.over.NoBirA, 
                    brain.acc$cellular_component, brain.acc$biological_process, brain.acc$molecular_function)
brn.nls_1 <- merge(test6, brn.nls, by = "Query")
write.csv(brn.nls_1, "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/NLS_predictions/brain_NLS_wGeneNames.csv")

# liver
test6 <- data.frame("Query" = liv.acc$id, "Gene" = liv.acc$geneSymbol, liv.acc$Log.P.Value.xBirA.over.NoBirA, liv.acc$adj.P.Val.xBirA.over.NoBirA, 
                    liv.acc$cellular_component, liv.acc$biological_process, liv.acc$molecular_function)
liv.nls_1 <- merge(test6, liv.nls, by = "Query")
write.csv(liv.nls_1, "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/NLS_predictions/liver_NLS_wGeneNames.csv")

# kidney
test6 <- data.frame("Query" = kid.acc$id, "Gene" = kid.acc$geneSymbol, kid.acc$Log.P.Value.xBirA.over.NoBirA, kid.acc$adj.P.Val.xBirA.over.NoBirA, 
                    kid.acc$cellular_component, kid.acc$biological_process, kid.acc$molecular_function)
kid.nls_1 <- merge(test6, kid.nls, by = "Query")
write.csv(kid.nls_1, "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/NLS_predictions/kidney_NLS_wGeneNames.csv")

# serum
test6 <- data.frame("Query" = serum_mlt.ft$id, "Gene" = serum_mlt.ft$geneSymbol, serum_mlt.ft$logFC.xBirA.Cre.over.NoBirA, serum_mlt.ft$adj.P.Val.xBirA.Cre.over.NoBirA, 
                    serum_mlt.ft$cellular_component, serum_mlt.ft$biological_process, serum_mlt.ft$molecular_function)
ser.nls_1 <- merge(test6, ser.nls, by = "Query")
write.csv(ser.nls_1, "/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/output/all_enriched/NLS_predictions/serum_NLS_wGeneNames.csv")


