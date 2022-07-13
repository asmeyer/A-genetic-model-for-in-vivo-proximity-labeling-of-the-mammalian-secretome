# 2021-07-13
# ASM 
# TR01: BirA/Sox2Cre paper 

# usage: find shared and unique protein lists from tissues to match venn diagram data for LC-MS/MS for enrichment score data for 
#tissues and log2FC for Serum. 

# last updated: 2021-07-11


################################################################## begin ####################################################################
library(ggvenn)

# read in filtered multimedian data; log2FC > 1.0 and adjusted p < 0.05
serum_mlt.ft <- read.csv("input/filtered_multmed/serum_filtered_logFC-adjP.csv") # not normalized, still need to filter 
serum_mlt.ft$geneSymbol <- with(serum_mlt.ft, ifelse(geneSymbol == "", accession_number, geneSymbol))


# read in ES data from RY 
brain_ES <- read.csv("input/ES/BRNdata.csv")
brain_ES$geneSymbol <- with(brain_ES, ifelse(geneSymbol == "", id, geneSymbol))
brain_ES.ft <- brain_ES[brain_ES$ES_score >= 5, ]

kid_ES <- read.csv("input/ES/KIDdata.csv")
kid_ES$geneSymbol <- with(kid_ES, ifelse(geneSymbol == "", id, geneSymbol))
kid_ES.ft <- kid_ES[kid_ES$ES_score >= 5, ]

liv_ES <- read.csv("input/ES/LIVdata.csv")
liv_ES$geneSymbol <- with(liv_ES, ifelse(geneSymbol == "", id, geneSymbol))
liv_ES.ft <- liv_ES[liv_ES$ES_score >= 5, ]



# all four tissues 
allT_wSer.ls.ES <- list("Brain" = brain_ES.ft$geneSymbol, "Liver" = liv_ES.ft$geneSymbol, 
                     "Kidney" = kid_ES.ft$geneSymbol, "Serum" = serum_mlt.ft$geneSymbol)

# no serum 
allTissues.ls.ES <- list("Brain" = brain_ES.ft$geneSymbol, "Liver" = liv_ES.ft$geneSymbol, "Kidney" = kid_ES.ft$geneSymbol)

# venn all tissues 
ggvenn(allTissues.ls.ES, fill_color = c("#257DCF", "#42853D", "#AA0A3C"), show_percentage = F)
ggvenn(allT_wSer.ls.ES, fill_color = c("#257DCF", "#42853D", "#AA0A3C", "grey30"), show_percentage = F)



# ASM function to get 4-tissue venn diagram intersects 
vennGroups <- function(tissue_list) {
  library(dplyr)
  library(ggplot2)
  library(ggvenn)
  A <- unlist(tissue_list[1])
  B <- unlist(tissue_list[2])
  C <- unlist(tissue_list[3])
  D <- unlist(tissue_list[4])
  allT <- intersect(intersect(intersect(A, B), C), D) # all shared 
  t1_2 <- setdiff(intersect(A, B), allT) %>% setdiff(intersect(intersect(A, B), C)) %>% setdiff(intersect(intersect(A, B), D)) # shared between 1 & 2 
  t1_3 <- setdiff(intersect(A, C), allT) %>% setdiff(intersect(intersect(A, B), C)) %>% setdiff(intersect(intersect(A, C), D))  # shared between 1 & 3
  t1_4 <- setdiff(intersect(A, D), allT) %>% setdiff(intersect(intersect(A, B), D)) %>% setdiff(intersect(intersect(A, C), D)) # shared between 1 & 4
  t2_3 <- setdiff(intersect(B, C), allT) %>% setdiff(intersect(intersect(A, B), C)) %>% setdiff(intersect(intersect(B, C), D)) # shared between 2 & 3 
  t2_4 <- setdiff(intersect(B, D), allT) %>% setdiff(intersect(intersect(A, B), D)) %>% setdiff(intersect(intersect(B, C), D)) # shared between 2 & 4 
  t3_4 <- setdiff(intersect(C, D), allT) %>% setdiff(intersect(intersect(A, C), D)) %>% setdiff(intersect(intersect(B, C), D)) # shared between 3 & 4
  t1_2_3 <- setdiff(intersect(intersect(A, B), C), allT) # shared between 1, 2, & 3 
  t1_2_4 <- setdiff(intersect(intersect(A, B), D), allT) # shared between 1, 2, & 4
  t1_3_4 <- setdiff(intersect(intersect(A, C), D), allT) # shared between 1, 3, & 4
  t2_3_4 <- setdiff(intersect(intersect(B, C), D), allT) # shared between 2, 3, & 4 
  t1 <- setdiff(A, B) %>% setdiff(C) %>% setdiff(D)
  t2 <- setdiff(B, A) %>% setdiff(C) %>% setdiff(D)
  t3 <- setdiff(C, A) %>% setdiff(B) %>% setdiff(D)
  t4 <- setdiff(D, A) %>% setdiff(B) %>% setdiff(C)
  vennOutput <- list("T1" = t1, "T2" = t2, "T3" = t3, "T4" = t4, "T1-2" = t1_2, "T1-3" = t1_3, "T1-4" = t1_4, "T2-3" = t2_3, "T2-4" = t2_4, 
                     "T3-4" = t3_4, "T1-2-3" = t1_2_3, "T1-2-4" = t1_2_4, "T1-3-4" = t1_3_4, "T2-3-4" = t2_3_4, "All" = allT)
  return(vennOutput)
  
}


# get venn intersects 
BX_allTs.ES.venn <- vennGroups(allT_wSer.ls.ES)

# rename cols 
BX_allTs.ES.venn <- data.frame(lapply(BX_allTs.ES.venn, "length<-", max(lengths(BX_allTs.ES.venn))))
colnames(BX_allTs.ES.venn)
colnames(BX_allTs.ES.venn) <- c("Brain", "Liver", "Kidney", "Serum", "Brain_Liver", "Brain_Kidney", "Brain_Serum", "Liver_Kidney", 
                             "Liver_Serum", "Kidney_Serum", "Brain_Liver_Kidney", "Brain_Liver_Serum", "Brain_Kidney_Serum", 
                             "Liver_Kidney_Serum", "All")
# save vennout lists 
write.csv(BX_allTs.ES.venn, "output/ES_lists/allTissue_venn_info_wSerum_ES.csv")
BX_allTs.ES.venn <- read.csv("output/ES_lists/allTissue_venn_info_wSerum_ES.csv")

# venn diagram 
plot1 <- ggvenn(allT_wSer.ls.ES, fill_color = c("#257DCF", "#42853D", "#AA0A3C", "grey30"), show_percentage = F, text_size = 10, set_name_size = 12)
ggsave(filename = "All_tissues_ES_venn_diagram.pdf", plot = plot1, path = "output/ES_lists/", width = 14, height = 12, units = "in", dpi = "retina")



## GO analysis 
lilibrary(clusterProfiler)
brain.ES.hits <- as.vector(BX_allTs.ES.venn$Brain[1:172])
brain.ES.GO <- brain_ES.ft[brain_ES.ft$geneSymbol %in% brain.ES.hits, ]
brain.ES.GO.SPTM <- brain.ES.GO[brain.ES.GO$signalP == TRUE | brain.ES.GO$TM, ]


## clusterProfiler
brain.ES.cc <- enrichGO(gene         = brain.ES.GO$id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'UNIPROT',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

brain.ES.bp <- enrichGO(gene         = brain.ES.GO$id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'UNIPROT',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

brain.ES.mf <- enrichGO(gene         = brain.ES.GO$id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'UNIPROT',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)


brain.ES.cc_plot <- dotplot(brain.ES.cc, font.size = 26)
brain.ES.bp_plot <- dotplot(brain.ES.bp, font.size = 26)
brain.ES.mf_plot <- dotplot(brain.ES.mf, font.size = 26)

ggsave(filename = "GO_brain-ES-cc.pdf", plot = brainT.cc_plot, path = "output/ES_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-ES-bp.pdf", plot = brainT.bp_plot, path = "output/ES_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-ES-mf.pdf", plot = brainT.mf_plot, path = "output/ES_lists/", width = 16, height = 12, units = "in", dpi = "retina")



brain.ES.SPTM.cc <- enrichGO(gene         = brain.ES.GO.SPTM$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brain.ES.SPTM.bp <- enrichGO(gene         = brain.ES.GO.SPTM$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brain.ES.SPTM.mf <- enrichGO(gene         = brain.ES.GO.SPTM$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)


brain.ES.SPTM.cc_plot <- dotplot(brain.ES.SPTM.cc, font.size = 26)
brain.ES.SPTM.bp_plot <- dotplot(brain.ES.SPTM.bp, font.size = 26)
brain.ES.SPTM.mf_plot <- dotplot(brain.ES.SPTM.mf, font.size = 26)

ggsave(filename = "GO_brain-ES-SPTM-cc.pdf", plot = brain.ES.SPTM.cc_plot, path = "output/ES_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-ES-SPTM-bp.pdf", plot = brain.ES.SPTM.bp_plot, path = "output/ES_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-ES-SPTM-mf.pdf", plot = brain.ES.SPTM.mf_plot, path = "output/ES_lists/", width = 16, height = 12, units = "in", dpi = "retina")


