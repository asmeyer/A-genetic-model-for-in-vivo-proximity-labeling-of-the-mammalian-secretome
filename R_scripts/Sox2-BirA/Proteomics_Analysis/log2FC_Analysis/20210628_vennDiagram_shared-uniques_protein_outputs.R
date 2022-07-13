# 2021-06-28
# ASM 
# TR01: BirA/Sox2Cre paper 

# usage: find shared and unique protein lists from tissues to match venn diagram data for LC-MS/MS. 

# last updated: 2021-10-17 (copied serum GO to new script for clarity)


################################################################## begin ####################################################################
library(ggvenn)

# read in filtered multimedian data; log2FC > 1.0 and adjusted p < 0.05
serum_mlt.ft <- read.csv("input/filtered_multmed/serum_filtered_logFC-adjP.csv") # not normalized, still need to filter 
serum_mlt.ft$geneSymbol <- with(serum_mlt.ft, ifelse(geneSymbol == "", accession_number, geneSymbol))



brain_mlt.ft <- read.csv("input/filtered_multmed/brain_filtered_logFC-adjP.csv")
brain_mlt.ft$geneSymbol <- with(brain_mlt.ft, ifelse(geneSymbol == "", accession_number, geneSymbol))

kid_mlt.ft <- read.csv("input/filtered_multmed/kid_filtered_logFC-adjP.csv")
kid_mlt.ft$geneSymbol <- with(kid_mlt.ft, ifelse(geneSymbol == "", accession_numbers, geneSymbol))

liv_mlt.ft <- read.csv("input/filtered_multmed/liv_filtered_logFC-adjP.csv")
liv_mlt.ft$geneSymbol <- with(liv_mlt.ft, ifelse(geneSymbol == "", id, geneSymbol))

allT_wSer.ls <- list("Brain" = brain_mlt.ft$geneSymbol, "Liver" = liv_mlt.ft$geneSymbol, 
                     "Kidney" = kid_mlt.ft$geneSymbol, "Serum" = serum_mlt.ft$geneSymbol)
allTissues.ls <- list("Brain" = brain_mlt.ft$geneSymbol, "Liver" = liv_mlt.ft$geneSymbol, "Kidney" = kid_mlt.ft$geneSymbol)

# venn all tissues 
ggvenn(allTissues.ls, fill_color = c("#257DCF", "#42853D", "#AA0A3C"), show_percentage = F)
ggvenn(allT_wSer.ls, fill_color = c("#257DCF", "#42853D", "#AA0A3C", "grey30"), show_percentage = F)


### shared all tissues 

### need to fill in blanks for geneSymbol in .flt files 


allTs.int <- intersect(intersect(brain_mlt.ft$geneSymbol, kid_mlt.ft$geneSymbol), liv_mlt.ft$geneSymbol) # shared all three; should be lenght 90
length(allTs.int)
View(allTs.int)
write.csv(allTs.int, "output/multmed_lists/Multmed_allTissues-intersect_noSerum.csv")

## brain_liv.int
brain_liv.int <- setdiff(intersect(brain_mlt.ft$geneSymbol, liv_mlt.ft$geneSymbol), allTs.int)
length(brain_liv.int)
write.csv(brain_liv.int, "output/multmed_lists/Multmed_brain-liv_intersect.csv")

## brain_kid.int
brain_kid.int <- setdiff(intersect(brain_mlt.ft$geneSymbol, kid_mlt.ft$geneSymbol), allTs.int)
length(brain_kid.int)
write.csv(brain_kid.int, "output/multmed_lists/Multmed_brain-kid_intersect.csv")

## kid_liv.int
kid_liv.int <- setdiff(intersect(kid_mlt.ft$geneSymbol, liv_mlt.ft$geneSymbol), allTs.int)
length(kid_liv.int)
write.csv(kid_liv.int, "output/multmed_lists/Multmed_kid-liv_intersect.csv")

## brain uniques 
brain.set <- setdiff(brain_mlt.ft$geneSymbol, kid_mlt.ft$geneSymbol) %>% 
  setdiff(liv_mlt.ft$geneSymbol)
length(brain.set)
write.csv(brain.set, "output/multmed_lists/Multmed_brain-setdiff.csv")

## liver uniques 
liv.set <- setdiff(liv_mlt.ft$geneSymbol, brain_mlt.ft$geneSymbol) %>% 
  setdiff(kid_mlt.ft$geneSymbol)
length(liv.set)
write.csv(liv.set, "output/multmed_lists/Multmed_liv-setdiff.csv")

## kidney uniques 
kid.set <- setdiff(kid_mlt.ft$geneSymbol, brain_mlt.ft$geneSymbol) %>% 
  setdiff(liv_mlt.ft$geneSymbol)
length(kid.set)
write.csv(kid.set, "output/multmed_lists/Multmed_kid-setdiff.csv")


## totals based on geneSymbol
# brain 51
# liv 43
# kid 381

## joints 
# brain-liv 2
# brain-kid 52
# liv-kid 48
# all 95


function(tissue_1, tissues_2, tissue_3, tissue_4) {
  
}


all.tissues.ls <- list("Brain" = brain_mlt.ft$geneSymbol, "Liver" = liv_mlt.ft$geneSymbol, 
                    "Kidney" = kid_mlt.ft$geneSymbol, "Serum" = serum_mlt.ft$geneSymbol)

str(tissue_list)

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







# function to get all unique and intersect values from venn diagram. Tissue order (in input list) directly corresponds to output comparisions.
# outputs are listed as T for tissue and numbers (1-4) denoting which tissue based on input list order. 
 
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


all.tissues.ls <- list("Brain" = brain_mlt.ft$geneSymbol, "Liver" = liv_mlt.ft$geneSymbol, 
                       "Kidney" = kid_mlt.ft$geneSymbol, "Serum" = serum_mlt.ft$geneSymbol)

BX_allTs.venn <- vennGroups(all.tissues.ls)
BX_allTs.venn <- data.frame(lapply(BX_allTs.venn, "length<-", max(lengths(BX_allTs.venn))))
colnames(BX_allTs.venn)
colnames(BX_allTs.venn) <- c("Brain", "Liver", "Kidney", "Serum", "Brain_Liver", "Brain_Kidney", "Brain_Serum", "Liver_Kidney", 
                             "Liver_Serum", "Kidney_Serum", "Brain_Liver_Kidney", "Brain_Liver_Serum", "Brain_Kidney_Serum", 
                             "Liver_Kidney_Serum", "All")
write.csv(BX_allTs.venn, "output/multmed_lists/allTissue_venn_info_wSerum.csv")
BX_allTs.venn <- read.csv("output/multmed_lists/allTissue_venn_info_wSerum.csv")
plot1 <- ggvenn(allT_wSer.ls, fill_color = c("#257DCF", "#42853D", "#AA0A3C", "grey30"), show_percentage = F, text_size = 10, set_name_size = 12)
ggsave(filename = "All_tissues_venn_diagram.pdf", plot = plot1, path = "output/multmed_lists/", width = 14, height = 12, units = "in", dpi = "retina")

str(allT_wSer.ls)
## GO analysis 
library(topGO)
library(rGREAT)

############################################## GO: log2FC & adj. p signs from each tissue ###################################################

library(clusterProfiler)


## clusterProfiler
brainT.cc <- enrichGO(gene         = brain_mlt.ft$accession_number,
                     OrgDb         = "org.Mm.eg.db",
                     keyType       = 'UNIPROT',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)

brainT.bp <- enrichGO(gene         = brain_mlt.ft$accession_number,
                     OrgDb         = "org.Mm.eg.db",
                     keyType       = 'UNIPROT',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)

brainT.mf <- enrichGO(gene         = brain_mlt.ft$accession_number,
                     OrgDb         = "org.Mm.eg.db",
                     keyType       = 'UNIPROT',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)


brainT.cc_plot <- dotplot(brainT.cc, font.size = 26)
brainT.bp_plot <- dotplot(brainT.bp, font.size = 26)
brainT.mf_plot <- dotplot(brainT.mf, font.size = 26)

ggsave(filename = "GO_brainT-cc.pdf", plot = brainT.cc_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brainT-bp.pdf", plot = brainT.bp_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brainT-mf.pdf", plot = brainT.mf_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")

livT.cc <- enrichGO(gene         = liv_mlt.ft$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

livT.bp <- enrichGO(gene         = liv_mlt.ft$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

livT.mf <- enrichGO(gene         = liv_mlt.ft$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)


livT.cc_plot <- dotplot(livT.cc, font.size = 26)
livT.bp_plot <- dotplot(livT.bp, font.size = 26)
livT.mf_plot <- dotplot(livT.mf, font.size = 26)

ggsave(filename = "GO_livT-cc.pdf", plot = livT.cc_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_livT-bp.pdf", plot = livT.bp_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_livT-mf.pdf", plot = livT.mf_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")


kidT.cc <- enrichGO(gene         = kid_mlt.ft$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

kidT.bp <- enrichGO(gene         = kid_mlt.ft$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

kidT.mf <- enrichGO(gene         = kid_mlt.ft$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)


kidT.cc_plot <- dotplot(kidT.cc, font.size = 26)
kidT.bp_plot <- dotplot(kidT.bp, font.size = 26)
kidT.mf_plot <- dotplot(kidT.mf, font.size = 26)

ggsave(filename = "GO_kidT-cc.pdf", plot = kidT.cc_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_kidT-bp.pdf", plot = kidT.bp_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_kidT-mf.pdf", plot = kidT.mf_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")


library(clusterProfiler)
serT.cc <- enrichGO(gene         = serum_mlt.ft$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

serT.bp <- enrichGO(gene         = serum_mlt.ft$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

serT.mf <- enrichGO(gene         = serum_mlt.ft$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)



serT.cc_plot <- dotplot(serT.cc, font.size = 26)
serT.bp_plot <- dotplot(serT.bp, font.size = 26)
serT.mf_plot <- dotplot(serT.mf, font.size = 26)

str(serT.cc_plot)
scale_y_continuous(limits=c(2,10))


ggsave(filename = "GO_serT-cc_v1.pdf", plot = serT.cc_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_serT-bp_v1.pdf", plot = serT.bp_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_serT-mf_v1.pdf", plot = serT.mf_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")





############################################## GO: unique to each tissue and shared from venn diagram ###################################################
shared.hits <- as.vector(BX_allTs.venn$All[1:37])
shared.GO<- brain_mlt.ft[brain_mlt.ft$geneSymbol %in% shared.hits, ]

brain.hits <- as.vector(BX_allTs.venn$Brain[1:41])
brain.GO <- brain_mlt.ft[brain_mlt.ft$geneSymbol %in% brain.hits, ]

liv.hits <- as.vector(BX_allTs.venn$Liver[1:26])
liv.GO <- liv_mlt.ft[liv_mlt.ft$geneSymbol %in% liv.hits, ]

kid.hits <- as.vector(BX_allTs.venn$Kidney[1:317])
kid.GO <- kid_mlt.ft[kid_mlt.ft$geneSymbol %in% kid.hits, ]

ser.hits <- as.vector(BX_allTs.venn$Serum[1:458])
ser.GO <- serum_mlt.ft[serum_mlt.ft$geneSymbol %in% ser.hits, ]

#### clusterProfiler
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

## clusterProfiler GO over representation analysis 
# cc = cellular component, bp = biological process, "MF" = molecular function
shared.cc <- enrichGO(gene         = shared.GO$accession_number,
                     OrgDb         = "org.Mm.eg.db",
                     keyType       = 'UNIPROT',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)

shared.bp <- enrichGO(gene         = shared.GO$accession_number,
                     OrgDb         = "org.Mm.eg.db",
                     keyType       = 'UNIPROT',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)

shared.mf <- enrichGO(gene         = shared.GO$accession_number,
                     OrgDb         = "org.Mm.eg.db",
                     keyType       = 'UNIPROT',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)


shared.cc_plot <- dotplot(shared.cc, font.size = 26)
shared.bp_plot <- dotplot(shared.bp, font.size = 26)
shared.mf_plot <- dotplot(shared.mf, font.size = 26)

ggsave(filename = "GO_shared-cc.pdf", plot = shared.cc_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_shared-bp.pdf", plot = shared.bp_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_shared-mf.pdf", plot = shared.mf_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")




brain.cc <- enrichGO(gene         = brain.GO$accession_number,
                 OrgDb         = "org.Mm.eg.db",
                 keyType       = 'UNIPROT',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

brain.bp <- enrichGO(gene         = brain.GO$accession_number,
                 OrgDb         = "org.Mm.eg.db",
                 keyType       = 'UNIPROT',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

brain.mf <- enrichGO(gene         = brain.GO$accession_number,
                     OrgDb         = "org.Mm.eg.db",
                     keyType       = 'UNIPROT',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)


brain.cc_plot <- dotplot(brain.cc, font.size = 26)
brain.bp_plot <- dotplot(brain.bp, font.size = 26)
brain.mf_plot <- dotplot(brain.mf, font.size = 26)

ggsave(filename = "GO_brain-cc.pdf", plot = brain.cc_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-bp.pdf", plot = brain.bp_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-mf.pdf", plot = brain.mf_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")

liv.cc <- enrichGO(gene         = liv.GO$id,
                     OrgDb         = "org.Mm.eg.db",
                     keyType       = 'UNIPROT',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)

liv.bp <- enrichGO(gene         = liv.GO$id,
                     OrgDb         = "org.Mm.eg.db",
                     keyType       = 'UNIPROT',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)

liv.mf <- enrichGO(gene         = liv.GO$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)


liv.cc_plot <- dotplot(liv.cc, font.size = 26)
liv.bp_plot <- dotplot(liv.bp, font.size = 26)
liv.mf_plot <- dotplot(liv.mf, font.size = 26)

ggsave(filename = "GO_liv-cc.pdf", plot = liv.cc_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_liv-bp.pdf", plot = liv.bp_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_liv-mf.pdf", plot = liv.mf_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")


kid.cc <- enrichGO(gene         = kid.GO$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

kid.bp <- enrichGO(gene         = kid.GO$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

kid.mf <- enrichGO(gene         = kid.GO$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)


kid.cc_plot <- dotplot(kid.cc, font.size = 26)
kid.bp_plot <- dotplot(kid.bp, font.size = 26)
kid.mf_plot <- dotplot(kid.mf, font.size = 26)

ggsave(filename = "GO_kid-cc.pdf", plot = kid.cc_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_kid-bp.pdf", plot = kid.bp_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_kid-mf.pdf", plot = kid.mf_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")



ser.cc <- enrichGO(gene         = ser.GO$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

ser.bp <- enrichGO(gene         = ser.GO$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

ser.mf <- enrichGO(gene         = ser.GO$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)



ser.cc_plot <- dotplot(ser.cc, font.size = 26)
ser.bp_plot <- dotplot(ser.bp, font.size = 26)
ser.mf_plot <- dotplot(ser.mf, font.size = 26)

ggsave(filename = "GO_ser-cc.pdf", plot = ser.cc_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_ser-bp.pdf", plot = ser.bp_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_ser-mf.pdf", plot = ser.mf_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")






####### permutation notes
# try combinations instead of permutations 


test <- list("A", "B", "C", "D")


library(combinat)
permn(letters[1:4])


library(gtools)
permutations(4, 4, letters[1:4])
permn(letters(1:i+1))
i = 0 
list[[i+1]]






### brain-kidney GO 
brk.hits <- as.vector(BX_allTs.venn$Brain_Kidney[1:41])
brk.GO <- brain_mlt.ft[brain_mlt.ft$geneSymbol %in% brk.hits, ]


brk.cc <- enrichGO(gene         = brk.GO$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

brk.bp <- enrichGO(gene         = brk.GO$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

brk.mf <- enrichGO(gene         = brk.GO$id,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)



(brk.cc_plot <- dotplot(brk.cc, font.size = 26))
(brk.bp_plot <- dotplot(brk.bp, font.size = 26))
(brk.mf_plot <- dotplot(brk.mf, font.size = 26))


