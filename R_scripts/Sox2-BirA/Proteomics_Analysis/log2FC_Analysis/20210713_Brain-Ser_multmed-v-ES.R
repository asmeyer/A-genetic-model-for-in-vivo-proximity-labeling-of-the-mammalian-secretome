# 2021-07-13
# ASM 
# TR01: BirA/Sox2Cre paper 

# usage: compare brain and serum unique proteins from multimedian data and ES data

# last updated: 2021-09-01


################################################################## begin ####################################################################
library(ggvenn)
library(ggplot2)

# read in filtered multimedian data; log2FC > 1.0 and adjusted p < 0.05
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


# ES data for each tissue 
brain.es <- read.csv("input/ES/BRN_unique.csv")
brain.es$geneSymbol <- with(brain.es, ifelse(geneSymbol == "", accession_number, geneSymbol))

serum.es <- read.csv("input/ES/Amanda_not_normalized_filtered_20201201.csv")
serum.es$Symbol_Original <- with(serum.es, ifelse(Symbol_Original  == "", accession, Symbol_Original))

liv.es <- read.csv("input/ES/LIV_unique.csv")
liv.es$Symbol_Original <- with(liv.es, ifelse(geneSymbol  == "", accession, geneSymbol))

kid.es <- read.csv("input/ES/KID_unique.csv")
kid.es$Symbol_Original <- with(kid.es, ifelse(geneSymbol  == "", accession, geneSymbol))

## all ES > 5
brain.es1 <- read.csv("input/ES/BRNdata.csv")
brain.es1$geneSymbol <- with(brain.es1, ifelse(geneSymbol == "", accession_number, geneSymbol))
brain.es1 <- brain.es1[brain.es1$ES_score >= 5, ]

liv.es1 <- read.csv("input/ES/LIVdata.csv")
liv.es1$geneSymbol <- with(liv.es1, ifelse(geneSymbol  == "", accession, geneSymbol))
liv.es1 <- liv.es1[liv.es1$ES_score >= 5, ]

kid.es1 <- read.csv("input/ES/KIDdata.csv")
kid.es1$geneSymbol <- with(kid.es1, ifelse(geneSymbol  == "", accession, geneSymbol))
kid.es1 <- kid.es1[kid.es1$ES_score >= 5, ]

# combine both for venn diagram using log2Fc & adjP < 0.05
brain.both1 <- list("Brain_MM" = brain_mlt.ft$geneSymbol, "Brain_ES" = brain.es1$geneSymbol)
serum.both1 <- list("Serum_MM" = serum_mlt.ft$geneSymbol, "Serum_ES" = serum.es$Symbol_Original)
liv.both1 <- list("Liver_MM" = liv_mlt.ft$geneSymbol, "Liver_ES" = liv.es1$geneSymbol)
kid.both1 <- list("Kidney_MM" = kid_mlt.ft$geneSymbol, "Kidney_ES" = kid.es1$geneSymbol)


# combine both for venn diagram using uniques 
brain.both <- list("Brain_MM" = brain.mm, "Brain_ES" = brain.es$geneSymbol)
serum.both <- list("Serum_MM" = serum.mm, "Serum_ES" = serum.es$Symbol_Original)
liv.both <- list("Liver_MM" = liv.mm, "Liver_ES" = liv.es$geneSymbol)
kid.both <- list("Kidney_MM" = kid.mm, "Kidney_ES" = kid.es$geneSymbol)

################################################################## MM vs. ES vennDiagrams ####################################################################
# brain
plot1 <- ggvenn(brain.both, fill_color = c("#257DCF", "#42853D"), show_percentage = F)
ggsave(filename = "Brain_multmed-ES_venn_diagram_unqiues.pdf", plot = plot1, path = "output/multmed_v_ES/", width = 14, height = 12, units = "in", dpi = "retina")
plot1 <- ggvenn(brain.both1, fill_color = c("#257DCF", "#42853D"), show_percentage = F)
ggsave(filename = "Brain_multmed-ES_venn_diagram.pdf", plot = plot1, path = "output/multmed_v_ES/", width = 14, height = 12, units = "in", dpi = "retina")

# serum
plot1 <- ggvenn(serum.both, fill_color = c("#257DCF", "#42853D"), show_percentage = F)
ggsave(filename = "serum_v2_multmed-ES_venn_diagram_uniques.pdf", plot = plot1, path = "output/multmed_v_ES/", width = 14, height = 12, units = "in", dpi = "retina")
plot1 <- ggvenn(serum.both1, fill_color = c("#257DCF", "#42853D"), show_percentage = F)
ggsave(filename = "serum_v2_multmed-ES_venn_diagram.pdf", plot = plot1, path = "output/multmed_v_ES/", width = 14, height = 12, units = "in", dpi = "retina")

# liver
plot1 <- ggvenn(liv.both, fill_color = c("#257DCF", "#42853D"), show_percentage = F)
ggsave(filename = "Liver_multmed-ES_venn_diagram_uniques.pdf", plot = plot1, path = "output/multmed_v_ES/", width = 14, height = 12, units = "in", dpi = "retina")
plot1 <- ggvenn(liv.both1, fill_color = c("#257DCF", "#42853D"), show_percentage = F)
ggsave(filename = "Liver_multmed-ES_venn_diagram.pdf", plot = plot1, path = "output/multmed_v_ES/", width = 14, height = 12, units = "in", dpi = "retina")

# kidney
plot1 <- ggvenn(kid.both, fill_color = c("#257DCF", "#42853D"), show_percentage = F)
ggsave(filename = "Kidney_multmed-ES_venn_diagram_uniques.pdf", plot = plot1, path = "output/multmed_v_ES/", width = 14, height = 12, units = "in", dpi = "retina")
plot1 <- ggvenn(kid.both1, fill_color = c("#257DCF", "#42853D"), show_percentage = F)
ggsave(filename = "Kidney_multmed-ES_venn_diagram.pdf", plot = plot1, path = "output/multmed_v_ES/", width = 14, height = 12, units = "in", dpi = "retina")


################################################################## Brain ####################################################################

# venn brain
plot1 <- ggvenn(brain.both, fill_color = c("#257DCF", "#42853D"), show_percentage = F)
ggsave(filename = "Brain_multmed-ES_venn_diagram.pdf", plot = plot1, path = "output/multmed_v_ES/", width = 14, height = 12, units = "in", dpi = "retina")

# shared and uniques 
brain.mm.uni1 <- setdiff(brain.mm, brain.es$geneSymbol)
brain.es.uni1 <- setdiff(brain.es$geneSymbol, brain.mm)
brain.shared1 <- intersect(brain.mm, brain.es$geneSymbol)

brain.mm.uni <- brain_mlt.ft[brain_mlt.ft$geneSymbol %in% brain.mm.uni1, ]
brain.es.uni <- brain.es[brain.es$geneSymbol %in% brain.es.uni1, ]
brain.shared <- brain_mlt.ft[brain_mlt.ft$geneSymbol %in% brain.shared1, ]

library(clusterProfiler)


## clusterProfiler
brain.mm.cc <- enrichGO(gene         = brain.mm.uni$id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'UNIPROT',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

brain.mm.bp <- enrichGO(gene         = brain.mm.uni$id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'UNIPROT',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

brain.mm.mf <- enrichGO(gene         = brain.mm.uni$id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'UNIPROT',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

brain.es.cc <- enrichGO(gene         = brain.es.uni$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brain.es.bp <- enrichGO(gene         = brain.es.uni$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brain.es.mf <- enrichGO(gene         = brain.es.uni$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brain.sh.cc <- enrichGO(gene         = brain.shared$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brain.sh.bp <- enrichGO(gene         = brain.shared$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brain.sh.mf <- enrichGO(gene         = brain.shared$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)


brain.mm.cc_plot <- dotplot(brain.mm.cc, font.size = 26)
brain.mm.bp_plot <- dotplot(brain.mm.bp, font.size = 26)
brain.mm.mf_plot <- dotplot(brain.mm.mf, font.size = 26)

brain.es.cc_plot <- dotplot(brain.es.cc, font.size = 26)
brain.es.bp_plot <- dotplot(brain.es.bp, font.size = 26)
brain.es.mf_plot <- dotplot(brain.es.mf, font.size = 26)

brain.sh.cc_plot <- dotplot(brain.sh.cc, font.size = 26)
brain.sh.bp_plot <- dotplot(brain.sh.bp, font.size = 26)
brain.sh.mf_plot <- dotplot(brain.sh.mf, font.size = 26)

ggsave(filename = "GO_brain-mm-cc.pdf", plot = brain.mm.cc_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-mm-bp.pdf", plot = brain.mm.bp_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-mm-mf.pdf", plot = brain.mm.mf_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")

ggsave(filename = "GO_brain-es-cc.pdf", plot = brain.es.cc_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-es-bp.pdf", plot = brain.es.bp_plot, path = "output/multmed_v_ES/", width = 24, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-es-mf.pdf", plot = brain.es.mf_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")

ggsave(filename = "GO_brain-sh-cc.pdf", plot = brain.sh.cc_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-sh-bp.pdf", plot = brain.sh.bp_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-sh-mf.pdf", plot = brain.sh.mf_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")

################################################################## Serum ####################################################################
# venn serum 
plot1 <- ggvenn(serum.both, fill_color = c("#257DCF", "#42853D"), show_percentage = F)
ggsave(filename = "Serum_multmed-ES_venn_diagram.pdf", plot = plot1, path = "output/multmed_v_ES/", width = 14, height = 12, units = "in", dpi = "retina")


# shared and uniques 
serum.mm.uni1 <- setdiff(serum.mm, serum.es$Symbol_Original)
serum.es.uni1 <- setdiff(serum.es$Symbol_Original, serum.mm)
serum.shared1 <- intersect(serum.mm, serum.es$Symbol_Original)

serum.mm.uni <- serum_mlt.ft[serum_mlt.ft$geneSymbol %in% serum.mm.uni1, ]
serum.es.uni <- serum.es[serum.es$Symbol_Original %in% serum.es.uni1, ]
serum.shared <- serum_mlt.ft[serum_mlt.ft$geneSymbol %in% serum.shared1, ]



## clusterProfiler
serum.mm.cc <- enrichGO(gene         = serum.mm.uni$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

serum.mm.bp <- enrichGO(gene         = serum.mm.uni$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

serum.mm.mf <- enrichGO(gene         = serum.mm.uni$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

serum.es.cc <- enrichGO(gene         = serum.es.uni$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

serum.es.bp <- enrichGO(gene         = serum.es.uni$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

serum.es.mf <- enrichGO(gene         = serum.es.uni$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

serum.sh.cc <- enrichGO(gene         = serum.shared$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

serum.sh.bp <- enrichGO(gene         = serum.shared$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

serum.sh.mf <- enrichGO(gene         = serum.shared$id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'UNIPROT',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)


serum.mm.cc_plot <- dotplot(serum.mm.cc, font.size = 26)
serum.mm.bp_plot <- dotplot(serum.mm.bp, font.size = 26)
serum.mm.mf_plot <- dotplot(serum.mm.mf, font.size = 26)

serum.es.cc_plot <- dotplot(serum.es.cc, font.size = 26)
serum.es.bp_plot <- dotplot(serum.es.bp, font.size = 26)
serum.es.mf_plot <- dotplot(serum.es.mf, font.size = 26)

serum.sh.cc_plot <- dotplot(serum.sh.cc, font.size = 26)
serum.sh.bp_plot <- dotplot(serum.sh.bp, font.size = 26)
serum.sh.mf_plot <- dotplot(serum.sh.mf, font.size = 26)

ggsave(filename = "GO_serum-mm-cc.pdf", plot = serum.mm.cc_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_serum-mm-bp.pdf", plot = serum.mm.bp_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_serum-mm-mf.pdf", plot = serum.mm.mf_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")

ggsave(filename = "GO_serum-es-cc.pdf", plot = serum.es.cc_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_serum-es-bp.pdf", plot = serum.es.bp_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_serum-es-mf.pdf", plot = serum.es.mf_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")

ggsave(filename = "GO_serum-sh-cc.pdf", plot = serum.sh.cc_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_serum-sh-bp.pdf", plot = serum.sh.bp_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_serum-sh-mf.pdf", plot = serum.sh.mf_plot, path = "output/multmed_v_ES/", width = 16, height = 12, units = "in", dpi = "retina")


# venn serum
ggvenn(allTissues.ls, fill_color = c("#257DCF", "#42853D", "#AA0A3C"), show_percentage = F)
ggvenn(allT_wSer.ls, fill_color = c("#257DCF", "#42853D", "#AA0A3C", "grey30"), show_percentage = F)
