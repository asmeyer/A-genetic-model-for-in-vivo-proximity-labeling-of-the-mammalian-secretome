# 2021-07-02
# ASM

# usage: clusterProfiler GO analysis example.



#### clusterProfiler
library(clusterProfiler)
library(org.Mm.eg.db) # for mouse GO analysis; change to org.Hs.eg.db for human
library(enrichplot)

## clusterProfiler GO over representation analysis 
# cc = cellular component, bp = biological process, "MF" = molecular function

# only need to change gene
# set to gene list of UNIPROT accession numbers 
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


(brain.cc_plot <- dotplot(brain.cc, font.size = 26))
(brain.bp_plot <- dotplot(brain.bp, font.size = 26))
(brain.mf_plot <- dotplot(brain.mf, font.size = 26))

ggsave(filename = "GO_brain-cc.pdf", plot = brain.cc_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-bp.pdf", plot = brain.bp_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brain-mf.pdf", plot = brain.mf_plot, path = "output/multmed_lists/", width = 16, height = 12, units = "in", dpi = "retina")














