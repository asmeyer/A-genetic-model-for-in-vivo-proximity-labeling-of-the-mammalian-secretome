# 20210922
# ASM

# usage: preliminary RNAseq analysis from WashU aligned data. Some parts from Louisa Gerhardt. 

# last updated: 2021-09-26

## Note: use v2 data. Re-did analysis to make Bir-sox2Cre/BirA (switched levels defaults). 

############################################################ begin #######################################################################
library(DESeq2)
library(ggplot2)

getwd()

# brain 
br92 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/92_br_1-25.AATAGGTAGG-TCAGCAAGCT/92_br_1-25.AATAGGTAGG-TCAGCAAGCT.gene_counts.txt", 
                   header = TRUE)
br93 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/93_br_1-25.CGATCCGCTT-GCAATGATGC/93_br_1-25.CGATCCGCTT-GCAATGATGC.gene_counts.txt", 
                   header = TRUE)
br94 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/94_br_1-25 .TCCGAGATGC-ATACGTGAGC/94_br_1-25 .TCCGAGATGC-ATACGTGAGC.gene_counts.txt", 
                 header = TRUE)
br95 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/95_br_1-25.GTCCACTCTC-ACGCCTCAAG/95_br_1-25.GTCCACTCTC-ACGCCTCAAG.gene_counts.txt", 
                 header = TRUE)

# liver
liv92 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/92_liv.AGGAGTCACT-TTCGACGCAG/92_liv.AGGAGTCACT-TTCGACGCAG.gene_counts.txt", 
                   header = TRUE)
liv93 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/93_liv.CAGTGGACCA-CTTCGGAATC/93_liv.CAGTGGACCA-CTTCGGAATC.gene_counts.txt", 
                   header = TRUE)
liv94 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/94_liv .ACTTAACGCG-GATGCGCATA/94_liv .ACTTAACGCG-GATGCGCATA.gene_counts.txt", 
                   header = TRUE)
liv95 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/95_liv.CCGCCAATAA-ACACATTGAC/95_liv.CCGCCAATAA-ACACATTGAC.gene_counts.txt", 
                   header = TRUE)

# kidney
kid92 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/92_kid.AACCATCAAG-CTCATGTCTA/92_kid.AACCATCAAG-CTCATGTCTA.gene_counts.txt", 
                   header = TRUE)
kid93 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/93_kid.AAGAACACTC-TATCCACAGC/93_kid.AAGAACACTC-TATCCACAGC.gene_counts.txt", 
                   header = TRUE)
kid94 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/94_kid.TCTTGCGTCT-ACACTCTCAA/94_kid.TCTTGCGTCT-ACACTCTCAA.gene_counts.txt", 
                   header = TRUE)
kid95 <- read.table("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/McMahon_s5467_MGI2053/95_kid.CACAGTAGGC-GAGTGTAATG/95_kid.CACAGTAGGC-GAGTGTAATG.gene_counts.txt", 
                   header = TRUE)

##################################################### DESeq all BirA/Sox2Cre to BirA #######################################################
# make count table 
allT.counts <- data.frame("Gene" = br92$Feature, "Br92" = br92$Count, "Br93" = br93$Count, "Br94" = br94$Count, "Br95" = br95$Count, 
                          "Liv92" = liv92$Count, "Liv93" = liv93$Count, "Liv94" = liv94$Count, "Liv95" = liv95$Count, 
                          "Kid92" = kid92$Count, "Kid93" = kid93$Count, "Kid94" = kid94$Count, "Kid95" = kid95$Count)

# make meta data table
allT.meta <- data.frame("id" = c("Br92", "Br93", "Br94", "Br95", "Liv92", "Liv93", "Liv94", "Liv95", "Kid92", "Kid93", "Kid94", "Kid95"), 
                        "genotype" = factor(rep(c("BirA", "BirA", "BirA-Sox2Cre", "BirA-Sox2Cre"), n = 4), levels = c("BirA-Sox2Cre", "BirA")), 
                        "tissue" = rep(c("brain", "liver", "kidney"), each = 4))
str(allT.meta$genotype)
str(allT.meta)

## deseq2 setup
allT.dds <- DESeqDataSetFromMatrix(countData = allT.counts, colData = allT.meta, design = ~genotype, tidy = TRUE)
str(allT.dds)

# differential expression analysis
allT.dds <- DESeq(allT.dds)
saveRDS(allT.dds, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/all_tissues_BirA_over_BirA-Sox2Cre_v2.Rds", compress = FALSE)
allT.dds <- readRDS("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/all_tissues_BirA_over_BirA-Sox2Cre_v2.Rds")

# results table
allT.res <- results(allT.dds)
head(results(allT.dds, tidy=TRUE))

# Write normalized gene counts to a .csv file
write.csv(x = as.data.frame(counts(All_datasets_DESeq, normalized = T)), 
          file = 'All_datasets_normalized_counts.csv')

# (Data has to be log transformed for that)
allT.dds_Rlog <- rlog(allT.dds, blind = FALSE)
all_pca_id_plot <- plotPCA(allT.dds_Rlog, intgroup = "id", ntop = nrow(allT.dds_Rlog)) + theme_classic()

# summary 
allSum <- summary(allT.mapped)
str(allSum)

# sorted by p-value
allT.res <- allT.res[order(allT.res$padj), ]
head(allT.res)

# add in geneSymbols
allT.mapped <- allT.res
allT.mapped$SYMBOL <- rownames(allT.mapped)
# use library AnnotationDB
allT.mapped$symbol <- mapIds(org.Mm.eg.db, keys = row.names(allT.res), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
allLog.mapped <- as.data.frame(allT.mapped@listData)

# filter based sig. p-value 
allsig <- allLog.mapped[allLog.mapped$padj <= 0.05, ]


summary(allT.res)
str(allT.res)


##### from Louisa 
#NA omit
allT_res.df <- na.omit(as.data.frame(allT.res))
range(allT_res.df$baseMean)


mcols(allT.dds, use.names=TRUE)

#add normalized count of each sample
norm_counts_df <- as.data.frame(counts(allT.dds, normalized = T))
norm_counts_df <- norm_counts_df[which(row.names(norm_counts_df) %in% row.names(allT_res.df)), ]

COL_data <- as.data.frame(allT.dds@colData)
brain <- COL_data[which(COL_data$tissue == "brain"), 4]
liver <- COL_data[which(COL_data$tissue == "liver"), 4]
kidney <- COL_data[which(COL_data$tissue == "kidney"), 4]

#Get mean for each condition
baseMean_brain <- rowMeans(norm_counts_df[ , 1:4])
baseMean_liver <- rowMeans(norm_counts_df[ , 5:8])
baseMean_kidney <- rowMeans(norm_counts_df[ , 9:12])

#combine df with normalized counts and mean for condition, sort df
norm_counts_df <- cbind(baseMean_brain, baseMean_liver, baseMean_kidney, norm_counts_df)

norm_counts_df <- norm_counts_df[ , c("baseMean_brain",as.character(brain),"baseMean_liver",as.character(liver), "baseMean_kidney",as.character(kidney))]

allT_res.df1 <- cbind(allT_res.df, norm_counts_df)
head(allT_res.df1)

### Annotating Gene Symbol/Description with BiomaRt
allT_res.df1 <- allT_res.df1 %>% tibble::rownames_to_column("ensembl_gene_id") 

library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

All_datasets.anno <- getBM(values = allT_res.df1$ensembl_gene_id, filters = "ensembl_gene_id", 
                           attributes = c("ensembl_gene_id", "entrezgene_id", "description", 'mgi_symbol'), mart = mart)

## use Louisa's nomenclature 
saveRDS(All_datasets.anno, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/All_datasets_anno_v2.Rds", compress = FALSE)
anno <- readRDS("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/All_datasets_anno_v2.Rds")

allT_res.df1.anno <- merge(anno, allT_res.df1, by="ensembl_gene_id")
head(allT_res.df1.anno)

#Keep only genes with absolute logFC >= 2 and p value adj < 0.05
df_DEG <- allT_res.df1.anno %>% arrange(-log2FoldChange) %>% filter(padj < 0.05 & abs(log2FoldChange) >= 2)

#reorder df_DEG and save csv
#df_DEG <- df_DEG[ , c(1,4,3,9,10,6,7,8,2,5,11:length(colnames(df_DEG)))]

write.csv(df_DEG, file = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/all_df_DEG_v2.csv", row.names = FALSE)

head(df_DEG, n = 20)

## Part 2-1: Volcano plot
## Create a column to indicate which genes to label
library(dplyr)
allT_res.df1.anno  <- mutate(allT_res.df1.anno , color = case_when(allT_res.df1.anno$log2FoldChange > 0 & allT_res.df1.anno$padj <0.05 ~ "Increased",
                                                             allT_res.df1.anno$log2FoldChange < 0 & allT_res.df1.anno$padj <0.05 ~ "Decreased",
                                                             allT_res.df1.anno$padj > 0.05 ~ "nonsignificant"))
allT_res.df1.anno <- allT_res.df1.anno %>% mutate(threshold = padj < 0.05) %>% arrange(padj) %>% mutate(volcanolabels = "")

#Top20
allT_res.df1.anno$volcanolabels[1:20] <- allT_res.df1.anno$mgi_symbol[1:20]

library(ggplot2)
library(ggrepel)
vol_p <-  ggplot(allT_res.df1.anno, aes(x = log2FoldChange, y = -log10(padj), color=color)) +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  geom_text_repel(aes(label = volcanolabels)) + scale_color_manual(name = "Directionality",
                                                                   values = c(Increased = "#D53E4F", Decreased = "#9ECAE1", nonsignificant = "grey50"))+
  xlab("log2(Fold Change)") + 
  ylab("-log10(Adj. p-value)") + theme_classic()
vol_p
ggsave(filename = "allTissues_DEG_volcano_plot_v2.pdf", plot = vol_p, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/", width = 6, height = 4, units = "in", dpi = "retina")
#"grey50", "#9ECAE1", "#C6DBEF", "#D53E4F"
## Part 2-2: Heatmap of top 15 up and top 15 downregulated genes
heatmap_genes <- df_DEG %>% arrange(padj) %>% top_n(n = -30, wt = padj)

# Gather 30 significant genes and make matrix
mat <- assay(allT.dds_Rlog)[heatmap_genes$ensembl_gene_id, ]
row.names(mat) <- heatmap_genes$mgi_symbol

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(allT.dds_Rlog)$Genotype), 
  row.names = colData(allT.dds_Rlog)$rownames
)

allT.dds_Rlog@colData$genotype

# Specify colors you want to annotate the columns by. Can add more groups like replicates here
ann_colors = list(Group = c(BirA = "grey50", BirA-Sox2Cre = "#9ECAE1"))

# Make Heatmap with pheatmap function.
## See more in documentation for customization
library(pheatmap)
plot_heatmap <- pheatmap(mat = mat, 
                         color = colorRampPalette(brewer.pal(9, "Blues"))(255), 
                         scale = "row", # Scale genes to Z-score (how many standard deviations)
                         annotation_colors = c("grey50", "#9ECAE1"),# Change the default colors of the annotations
                         fontsize = 6.5, # Make fonts smaller
                         cellwidth = 15, # Make the cells wider
                         show_colnames = T)



##Single gene plot
# Get gene with highest expression
top_gene <- allT_res.df1.anno$ensembl_gene_id[which.max(allT_res.df1.anno$log2FoldChange)]
plot_title <- allT_res.df1.anno$mgi_symbol[which(allT_res.df1.anno$ensembl_gene_id == top_gene)]

# Plot single gene
plotCounts(dds = allT.dds, 
           gene = top_gene, 
           intgroup = "genotype", 
           normalized = T, 
           transform = T, main =plot_title)



############################################################## all tissues PCA #############################################################
# (Data has to be log transformed for that)
allT.dds_Rlog <- rlog(allT.dds, blind = FALSE)
all_pca_id_plot <- plotPCA(allT.dds_Rlog, intgroup = "id", ntop = nrow(allT.dds_Rlog)) + theme_classic()

#### Genotype 
# get PCA data out 
library(RColorBrewer)
plot <- plotPCA(allT.dds_Rlog, intgroup = "genotype", ntop = nrow(allT.dds_Rlog)) + theme_classic() 
str(plot)
plot1 <- plot$data
plot1$tissue <- rep(c("Brain", "Liver", "Kidney"), each = 4)

# get PCA percentages 
percentVar1 <- plot$labels

all_pca_gentoype_plot <- ggplot(plot1, aes(x= PC1, y = PC2)) +
  geom_point(size = 7, alpha = 0.8, aes(color = genotype, shape = tissue)) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c("grey70", "#9ECAE1")) +
  xlab(paste0(percentVar1[2])) +
  ylab(paste0(percentVar1[1])) +
  coord_fixed() +
  labs(color = "Genotype", shape = "Tissue") +
  ggtitle("PCA with Scaled data") + theme(plot.title = element_text(hjust = 0.5))

#### Tissue
# get PCA data out 
library(RColorBrewer)
plot <- plotPCA(allT.dds_Rlog, intgroup = "tissue", ntop = nrow(allT.dds_Rlog)) + theme_classic() 
str(plot$data)
plot1 <- plot$data
plot1$Genotype <- rep(c("BirA", "BirA", "BirA-Sox2Cre", "BirA-Sox2Cre"), n = 4)

# get PCA percentages 
percentVar1 <- plot$labels

all_pca_tissue_plot <- ggplot(plot1, aes(x= PC1, y = PC2)) +
  geom_point(size = 7, alpha = 0.8, aes(color = tissue, shape = Genotype)) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c("#868686FF", "#628BB3", "#F6E0A3")) +
  xlab(paste0(percentVar1[2])) +
  ylab(paste0(percentVar1[1])) +
  coord_fixed() +
  labs(color = "Tissue", shape = "Genotype") +
  ggtitle("PCA with Scaled data") + theme(plot.title = element_text(hjust = 0.5))


### save 
ggsave(filename = "All_tissues_PCA_by_Genotype_v2.pdf", plot = all_pca_gentoype_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/", width = 14, height = 12, units = "in", dpi = "retina")
ggsave(filename = "All_tissues_PCA_by_Tissue_v2.pdf", plot = all_pca_tissue_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/",  width = 14, height = 12, units = "in", dpi = "retina")





########################################################## GO analysis all tissues DEGs ##################################################################

## clusterProfiler
# gene ontology for brain BirA enriched genes over BirA/Sox2Cre 
library(clusterProfiler)
library(ggplot2)
df_DEG <- read.csv("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/all_df_DEG_v2.csv")
all.res.df1.up_regulated <- df_DEG[which(df_DEG$log2FoldChange > 0), ]


allRNA.cc <- enrichGO(gene         = all.res.df1.up_regulated$ensembl_gene_id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

allRNA.bp <- enrichGO(gene         = all.res.df1.up_regulated$ensembl_gene_id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

allRNA.mf <- enrichGO(gene         = all.res.df1.up_regulated$ensembl_gene_id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)


# no enriched terms found 
allRNA.cc_plot <- dotplot(allRNA.cc, font.size = 16)
allRNA.bp_plot <- dotplot(allRNA.bp, font.size = 16)
allRNA.mf_plot <- dotplot(allRNA.mf, font.size = 16)

ggsave(filename = "GO_allRNA_pos_LG_cc_plot_v2.pdf", plot = allRNA.cc_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_allRNA_pos_LG_bp_plot_v2.pdf", plot = allRNA.bp_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_allRNA_pos_LG_mf_plot_v2.pdf", plot = allRNA.mf_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/",  width = 16, height = 12, units = "in", dpi = "retina")



### negative
# BirA/Sox2Cre enriched genes over BirA
df_DEG <- read.csv("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps/all_df_DEG_v2.csv")
all.res.df1.down_regulated <- df_DEG[which(df_DEG$log2FoldChange < 0), ]


allRNA.cc <- enrichGO(gene         = all.res.df1.down_regulated$ensembl_gene_id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

allRNA.bp <- enrichGO(gene         = all.res.df1.down_regulated$ensembl_gene_id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

allRNA.mf <- enrichGO(gene         = all.res.df1.down_regulated$ensembl_gene_id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)


# no enriched terms found 
allRNA.cc_plot <- dotplot(allRNA.cc, font.size = 16)
allRNA.bp_plot <- dotplot(allRNA.bp, font.size = 16)
allRNA.mf_plot <- dotplot(allRNA.mf, font.size = 16)

ggsave(filename = "GO_allRNA_neg_LG_cc_plot_v2.pdf", plot = allRNA.cc_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps//", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_allRNA_neg_LG_bp_plot_v2.pdf", plot = allRNA.bp_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps//", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_allRNA_neg_LG_mf_plot_v2.pdf", plot = allRNA.mf_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/all_tissues_comps//",  width = 16, height = 12, units = "in", dpi = "retina")




##################################################### DESeq Brain BirA/Sox2Cre to BirA #######################################################
# make count table 
brain.counts <- data.frame("Gene" = br92$Feature, "Br92" = br92$Count, "Br93" = br93$Count, "Br94" = br94$Count, "Br95" = br95$Count)

# make meta data table
brain.meta <- data.frame("id" = c("Br92", "Br93", "Br94", "Br95"), 
                        "genotype" = factor(c("BirA", "BirA", "BirA-Sox2Cre", "BirA-Sox2Cre"), levels = c("BirA-Sox2Cre", "BirA")), 
                        "tissue" = rep(c("brain"), n = 4))



## deseq2 setup
brain.dds <- DESeqDataSetFromMatrix(countData = brain.counts, colData = brain.meta, design = ~genotype, tidy = TRUE)

# differential expression analysis
brain.dds <- DESeq(brain.dds)
saveRDS(brain.dds, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/brain_BirA_over_BirA-Sox2Cre.Rds", compress = FALSE)
brain.dds <- readRDS("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/brain_BirA_over_BirA-Sox2Cre.Rds")


### PCA 
# (Data has to be log transformed for that)
brain.dds_Rlog <- rlog(brain.dds, blind = FALSE)

plotPCA(brain.dds_Rlog, intgroup = "id", ntop = nrow(brain.dds_Rlog)) + theme_classic()
plotPCA(brain.dds_Rlog, intgroup = "genotype", ntop = nrow(brain.dds_Rlog)) + theme_classic()


# results table
brain.res <- results(brain.dds)
head(results(brain.dds, tidy=TRUE))

# summary 
summary(brain.res)


# sorted by p-value
brain.res <- brain.res[order(brain.res$padj), ]
head(brain.res)

# add in geneSymbols
brain.mapped <- brain.res
brain.mapped$SYMBOL <- rownames(brain.mapped)
# use library AnnotationDB
brain.mapped$symbol <- mapIds(org.Mm.eg.db, keys = row.names(brain.res), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
brainLog.mapped <- as.data.frame(brain.mapped@listData)

# filter based sig. p-value 
brainsig <- brainLog.mapped[brainLog.mapped$padj <= 0.05, ]
write.csv(brainsig, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/brain_sig_BirA-BirA-Sox2Cre_v2.csv")

brainsig <- read.csv("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/brain_sig_BirA-BirA-Sox2Cre_v2.csv")
brain.acc <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/input/TR01_Brain_MedianNormalized_Multimedian.csv")

# BirA enriched genes over BirA/Sox2Cre
brainsig_pos <- brainsig[brainsig$log2FoldChange > 0.00, ]

##### from Louisa 
#NA omit
brain.res.df <- na.omit(as.data.frame(brain.res))
range(brain.res.df$baseMean)


mcols(brain.dds, use.names=TRUE)

#add normalized count of each sample
norm_counts_df <- as.data.frame(counts(brain.dds, normalized = T))
norm_counts_df <- norm_counts_df[which(row.names(norm_counts_df) %in% row.names(brain.res.df)), ]

COL_data <- as.data.frame(brain.dds@colData)
BirA <- COL_data[which(COL_data$tissue == "BirA"), 4]
BirA_Sox2Cre <- COL_data[which(COL_data$tissue == "BirA-Sox2Cre"), 4]


#Get mean for each condition
baseMean_BirA <- rowMeans(norm_counts_df[ , 1:2])
baseMean_BirA_Sox2Cre <- rowMeans(norm_counts_df[ , 3:4])


#combine df with normalized counts and mean for condition, sort df
norm_counts_df <- cbind(baseMean_BirA, baseMean_BirA_Sox2Cre, norm_counts_df)


brain.res.df1 <- cbind(brain.res.df, norm_counts_df)
head(brain.res.df1)

### Annotating Gene Symbol/Description with BiomaRt
brain.res.df1 <- brain.res.df1 %>% tibble::rownames_to_column("ensembl_gene_id") 

library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

brain_datasets.anno <- getBM(values = brain.res.df1$ensembl_gene_id, filters = "ensembl_gene_id", 
                           attributes = c("ensembl_gene_id", "entrezgene_id", "description", 'mgi_symbol'), mart = mart)

## use Louisa's nomenclature 
saveRDS(brain_datasets.anno, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/brain_datasets_anno_v2.Rds", compress = FALSE)
anno <- readRDS("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/brain_datasets_anno_v2.Rds")

brain.res.df1.anno <- merge(anno, brain.res.df1, by="ensembl_gene_id")
head(brain.res.df1.anno)

#Keep only genes with absolute logFC >= 2 and p value adj < 0.05
df_DEG <- brain.res.df1.anno %>% arrange(-log2FoldChange) %>% filter(padj < 0.05 & abs(log2FoldChange) >= 2)

#reorder df_DEG and save csv
#df_DEG <- df_DEG[ , c(1,4,3,9,10,6,7,8,2,5,11:length(colnames(df_DEG)))]

write.csv(df_DEG, file = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/Brain_df_DEG_v2.csv", row.names = FALSE)

head(df_DEG, n = 20)

## Part 2-1: Volcano plot
## Create a column to indicate which genes to label
library(dplyr)
brain.res.df1.anno  <- mutate(brain.res.df1.anno , color = case_when(brain.res.df1.anno$log2FoldChange > 0 & brain.res.df1.anno$padj <0.05 ~ "Increased",
                                                                    brain.res.df1.anno$log2FoldChange < 0 & brain.res.df1.anno$padj <0.05 ~ "Decreased",
                                                                    brain.res.df1.anno$padj > 0.05 ~ "nonsignificant"))
brain.res.df1.anno <- brain.res.df1.anno %>% mutate(threshold = padj < 0.05) %>% arrange(padj) %>% mutate(volcanolabels = "")

#Top20
brain.res.df1.anno$volcanolabels[1:20] <- brain.res.df1.anno$mgi_symbol[1:20]

library(ggplot2)
library(ggrepel)
vol_p <-  ggplot(brain.res.df1.anno, aes(x = log2FoldChange, y = -log10(padj), color=color)) +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  geom_text_repel(aes(label = volcanolabels)) + scale_color_manual(name = "Directionality",
                                                                   values = c(Increased = "#D53E4F", Decreased = "#9ECAE1", nonsignificant = "grey50"))+
  xlab("log2(Fold Change)") + 
  ylab("-log10(Adj. p-value)") + theme_classic()
vol_p
ggsave(filename = "brain_DEG_volcano_plot_v2.pdf", plot = vol_p, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain", width = 6, height = 4, units = "in", dpi = "retina")

#"grey50", "#9ECAE1", "#C6DBEF", "#D53E4F"
## Part 2-2: Heatmap of top 15 up and top 15 downregulated genes
heatmap_genes <- df_DEG %>% arrange(padj) %>% top_n(n = -30, wt = padj)

# Gather 30 significant genes and make matrix
mat <- assay(brain.dds_Rlog)[heatmap_genes$ensembl_gene_id, ]
row.names(mat) <- heatmap_genes$mgi_symbol

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(brain.dds_Rlog)$Genotype), 
  row.names = colData(brain.dds_Rlog)$rownames
)

brain.dds_Rlog@colData$genotype

# Specify colors you want to annotate the columns by. Can add more groups like replicates here
ann_colors = list(Group = c(BirA = "grey50", BirA-Sox2Cre = "#9ECAE1"))

# Make Heatmap with pheatmap function.
## See more in documentation for customization
library(pheatmap)
plot_heatmap <- pheatmap(mat = mat, 
                         color = colorRampPalette(brewer.pal(9, "Blues"))(255), 
                         scale = "row", # Scale genes to Z-score (how many standard deviations)
                         annotation_colors = c("grey50", "#9ECAE1"),# Change the default colors of the annotations
                         fontsize = 6.5, # Make fonts smaller
                         cellwidth = 15, # Make the cells wider
                         show_colnames = T)



##Single gene plot
# Get gene with highest expression
top_gene <- brain.res.df1.anno$ensembl_gene_id[which.max(brain.res.df1.anno$log2FoldChange)]
plot_title <- brain.res.df1.anno$mgi_symbol[which(brain.res.df1.anno$ensembl_gene_id == top_gene)]

# Plot single gene
plotCounts(dds = brain.dds, 
           gene = top_gene, 
           intgroup = "genotype", 
           normalized = T, 
           transform = T, main =plot_title)









############################################################## GO analysis Brain ######################################################################
## clusterProfiler
# gene ontology for brain BirA enriched genes over BirA/Sox2Cre 
library(clusterProfiler)
library(ggplot2)
df_DEG <- read.csv("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/brain_df_DEG_v2.csv")
brain.res.df1.up_regulated <- df_DEG[which(df_DEG$log2FoldChange > 0), ]

brainRNA.cc <- enrichGO(gene         = brainsig_pos$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

brainRNA.bp <- enrichGO(gene         = brainsig_pos$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

brainRNA.mf <- enrichGO(gene         = brainsig_pos$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

brainRNA.cc <- enrichGO(gene         = brain.res.df1.up_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

brainRNA.bp <- enrichGO(gene         = brain.res.df1.up_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

brainRNA.mf <- enrichGO(gene         = brain.res.df1.up_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)


# no enriched terms found 
brainRNA.cc_plot <- dotplot(brainRNA.cc, font.size = 16)
brainRNA.bp_plot <- dotplot(brainRNA.bp, font.size = 16)
brainRNA.mf_plot <- dotplot(brainRNA.mf, font.size = 16)

ggsave(filename = "GO_brainRNA_pos_LG_cc_plot_v2.pdf", plot = brainRNA.cc_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brainRNA_pos_LG_bp_plot_v2.pdf", plot = brainRNA.bp_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brainRNA_pos_LG_mf_plot_V2.pdf", plot = brainRNA.mf_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/",  width = 16, height = 12, units = "in", dpi = "retina")



### negative
# BirA/Sox2Cre enriched genes over BirA
brainsig_neg <- brainsig[brainsig$log2FoldChange < 0.00, ]
df_DEG <- read.csv("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/brain_df_DEG_v2.csv")
brain.res.df1.down_regulated <- df_DEG[which(df_DEG$log2FoldChange < 0), ]



## clusterProfiler
# gene ontology for brain BirA/Sox2Cre enriched genes over BirA 
library(clusterProfiler)
library(ggplot2)

brainRNA.cc <- enrichGO(gene         = brainsig_neg$SYMBOL,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brainRNA.bp <- enrichGO(gene         = brainsig_neg$SYMBOL,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brainRNA.mf <- enrichGO(gene         = brainsig_neg$SYMBOL,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brainRNA.cc <- enrichGO(gene         = brain.res.df1.down_regulated$ensembl_gene_id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brainRNA.bp <- enrichGO(gene         = brain.res.df1.down_regulated$ensembl_gene_id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

brainRNA.mf <- enrichGO(gene         = brain.res.df1.down_regulated$ensembl_gene_id,
                        OrgDb         = "org.Mm.eg.db",
                        keyType       = 'ENSEMBL',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)


# no enriched terms found 
brainRNA.cc_plot <- dotplot(brainRNA.cc, font.size = 16)
brainRNA.bp_plot <- dotplot(brainRNA.bp, font.size = 16)
brainRNA.mf_plot <- dotplot(brainRNA.mf, font.size = 16)

ggsave(filename = "GO_brainRNA_neg_LG_cc_plot_v2.pdf", plot = brainRNA.cc_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brainRNA_neg_LG_bp_plot_v2.pdf", plot = brainRNA.bp_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_brainRNA_neg_LG_mf_plot_v2.pdf", plot = brainRNA.mf_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/brain/",  width = 16, height = 12, units = "in", dpi = "retina")


##################################################### DESeq Liver BirA/Sox2Cre to BirA #######################################################
# make count table 
liv.counts <- data.frame("Gene" = liv92$Feature, "Liv92" = liv92$Count, "Liv93" = liv93$Count, "Liv94" = liv94$Count, "Liv95" = liv95$Count)

# make meta data table
liv.meta <- data.frame("id" = c("Liv92", "Liv93", "Liv94", "Liv95"), 
                       "genotype" = factor(c("BirA", "BirA", "BirA-Sox2Cre", "BirA-Sox2Cre"), levels = c("BirA-Sox2Cre", "BirA")), 
"tissue" = rep(c("liver"), n = 4))

## deseq2 setup
liv.dds <- DESeqDataSetFromMatrix(countData = liv.counts, colData = liv.meta, design = ~genotype, tidy = TRUE)

# differential expression analysis
liv.dds <- DESeq(liv.dds)
saveRDS(liv.dds, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/liver_BirA_over_BirA-Sox2Cre_v2.Rds", compress = FALSE)
liv.dds <- readRDS("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/liver_BirA_over_BirA-Sox2Cre_v2.Rds")

# (Data has to be log transformed for that)
liv.dds_Rlog <- rlog(liv.dds, blind = FALSE)

plotPCA(liv.dds_Rlog, intgroup = "id", ntop = nrow(liv.dds_Rlog)) + theme_classic()
plotPCA(liv.dds_Rlog, intgroup = "genotype", ntop = nrow(liv.dds_Rlog)) + theme_classic()


# results table
liv.res <- results(liv.dds)
head(results(liv.dds, tidy=TRUE))

# summary 
summary(liv.res)


# sorted by p-value
liv.res <- liv.res[order(liv.res$padj), ]
head(liv.res)

# add in geneSymbols
liv.mapped <- liv.res
liv.mapped$SYMBOL <- rownames(liv.mapped)
# use library AnnotationDB
liv.mapped$symbol <- mapIds(org.Mm.eg.db, keys = row.names(liv.res), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
livLog.mapped <- as.data.frame(liv.mapped@listData)

# filter based sig. p-value 
livsig <- livLog.mapped[livLog.mapped$padj <= 0.05, ]
write.csv(livsig, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/liver_sig_BirA-BirA-Sox2Cre_v2.csv")


livsig <- read.csv("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/liver_sig_BirA-BirA-Sox2Cre_v2.csv")
liv.acc <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/input/TR01_Liver_MedianNormalized_Multimedian.csv")

##### from Louisa 
#NA omit
liv.res.df <- na.omit(as.data.frame(liv.res))
range(liv.res.df$baseMean)


mcols(liv.dds, use.names=TRUE)

#add normalized count of each sample
norm_counts_df <- as.data.frame(counts(liv.dds, normalized = T))
norm_counts_df <- norm_counts_df[which(row.names(norm_counts_df) %in% row.names(liv.res.df)), ]

COL_data <- as.data.frame(liv.dds@colData)
BirA <- COL_data[which(COL_data$tissue == "BirA"), 4]
BirA_Sox2Cre <- COL_data[which(COL_data$tissue == "BirA-Sox2Cre"), 4]


#Get mean for each condition
baseMean_BirA <- rowMeans(norm_counts_df[ , 1:2])
baseMean_BirA_Sox2Cre <- rowMeans(norm_counts_df[ , 3:4])


#combine df with normalized counts and mean for condition, sort df
norm_counts_df <- cbind(baseMean_BirA, baseMean_BirA_Sox2Cre, norm_counts_df)


liv.res.df1 <- cbind(liv.res.df, norm_counts_df)
head(liv.res.df1)

### Annotating Gene Symbol/Description with BiomaRt
liv.res.df1 <- liv.res.df1 %>% tibble::rownames_to_column("ensembl_gene_id") 

library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

liv_datasets.anno <- getBM(values = liv.res.df1$ensembl_gene_id, filters = "ensembl_gene_id", 
                             attributes = c("ensembl_gene_id", "entrezgene_id", "description", 'mgi_symbol'), mart = mart)

## use Louisa's nomenclature 
saveRDS(liv_datasets.anno, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/liv_datasets_anno_v2.Rds", compress = FALSE)
anno <- readRDS("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/liv_datasets_anno_v2.Rds")

liv.res.df1.anno <- merge(anno, liv.res.df1, by="ensembl_gene_id")
head(liv.res.df1.anno)

#Keep only genes with absolute logFC >= 2 and p value adj < 0.05
df_DEG <- liv.res.df1.anno %>% arrange(-log2FoldChange) %>% filter(padj < 0.05 & abs(log2FoldChange) >= 2)

#reorder df_DEG and save csv
#df_DEG <- df_DEG[ , c(1,4,3,9,10,6,7,8,2,5,11:length(colnames(df_DEG)))]

write.csv(df_DEG, file = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/liv_df_DEG_v2.csv", row.names = FALSE)

head(df_DEG, n = 20)

## Part 2-1: Volcano plot
## Create a column to indicate which genes to label
library(dplyr)
liv.res.df1.anno  <- mutate(liv.res.df1.anno , color = case_when(liv.res.df1.anno$log2FoldChange > 0 & liv.res.df1.anno$padj <0.05 ~ "Increased",
                                                                     liv.res.df1.anno$log2FoldChange < 0 & liv.res.df1.anno$padj <0.05 ~ "Decreased",
                                                                     liv.res.df1.anno$padj > 0.05 ~ "nonsignificant"))
liv.res.df1.anno <- liv.res.df1.anno %>% mutate(threshold = padj < 0.05) %>% arrange(padj) %>% mutate(volcanolabels = "")

#Top20
liv.res.df1.anno$volcanolabels[1:20] <- liv.res.df1.anno$mgi_symbol[1:20]

library(ggplot2)
library(ggrepel)
vol_p <-  ggplot(liv.res.df1.anno, aes(x = log2FoldChange, y = -log10(padj), color=color)) +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  geom_text_repel(aes(label = volcanolabels)) + scale_color_manual(name = "Directionality",
                                                                   values = c(Increased = "#D53E4F", Decreased = "#9ECAE1", nonsignificant = "grey50"))+
  xlab("log2(Fold Change)") + 
  ylab("-log10(Adj. p-value)") + theme_classic()
vol_p
ggsave(filename = "liv_DEG_volcano_plot_v2.pdf", plot = vol_p, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver", width = 6, height = 4, units = "in", dpi = "retina")

#"grey50", "#9ECAE1", "#C6DBEF", "#D53E4F"
## Part 2-2: Heatmap of top 15 up and top 15 downregulated genes
heatmap_genes <- df_DEG %>% arrange(padj) %>% top_n(n = -30, wt = padj)

# Gather 30 significant genes and make matrix
mat <- assay(liv.dds_Rlog)[heatmap_genes$ensembl_gene_id, ]
row.names(mat) <- heatmap_genes$mgi_symbol

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(liv.dds_Rlog)$Genotype), 
  row.names = colData(liv.dds_Rlog)$rownames
)

liv.dds_Rlog@colData$genotype

# Specify colors you want to annotate the columns by. Can add more groups like replicates here
ann_colors = list(Group = c(BirA = "grey50", BirA-Sox2Cre = "#9ECAE1"))

# Make Heatmap with pheatmap function.
## See more in documentation for customization
library(pheatmap)
plot_heatmap <- pheatmap(mat = mat, 
                         color = colorRampPalette(brewer.pal(9, "Blues"))(255), 
                         scale = "row", # Scale genes to Z-score (how many standard deviations)
                         annotation_colors = c("grey50", "#9ECAE1"),# Change the default colors of the annotations
                         fontsize = 6.5, # Make fonts smaller
                         cellwidth = 15, # Make the cells wider
                         show_colnames = T)



##Single gene plot
# Get gene with highest expression
top_gene <- liv.res.df1.anno$ensembl_gene_id[which.max(liv.res.df1.anno$log2FoldChange)]
plot_title <- liv.res.df1.anno$mgi_symbol[which(liv.res.df1.anno$ensembl_gene_id == top_gene)]

# Plot single gene
plotCounts(dds = liv.dds, 
           gene = top_gene, 
           intgroup = "genotype", 
           normalized = T, 
           transform = T, main =plot_title)



################################################################ GO analysis liver ###########################################################################

# BirA enriched genes over BirA/Sox2Cre
livsig_pos <- livsig[livsig$log2FoldChange > 0.00, ]
df_DEG <- read.csv("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/liv_df_DEG_v2.csv")
liv.res.df1.up_regulated <- df_DEG[which(df_DEG$log2FoldChange > 0), ]


## clusterProfiler
# gene ontology for liver BirA enriched genes over BirA/Sox2Cre 
library(clusterProfiler)
library(ggplot2)

livRNA.cc <- enrichGO(gene         = livsig_pos$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

livRNA.bp <- enrichGO(gene         = livsig_pos$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

livRNA.mf <- enrichGO(gene         = livsig_pos$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

livRNA.cc <- enrichGO(gene         = liv.res.df1.up_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

livRNA.bp <- enrichGO(gene         = liv.res.df1.up_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

livRNA.mf <- enrichGO(gene         = liv.res.df1.up_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)


# no enriched terms found 
livRNA.cc_plot <- dotplot(livRNA.cc, font.size = 16)
livRNA.bp_plot <- dotplot(livRNA.bp, font.size = 16)
livRNA.mf_plot <- dotplot(livRNA.mf, font.size = 16)

ggsave(filename = "GO_livRNA_pos_LG_cc_plot_v2.pdf", plot = livRNA.cc_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_livRNA_pos_LG_bp_plot_v2.pdf", plot = livRNA.bp_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_livRNA_pos_LG_mf_plot_v2.pdf", plot = livRNA.mf_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/",  width = 16, height = 12, units = "in", dpi = "retina")



#### negative 
# BirA/Sox2Cre enriched genes over BirA
livsig_neg <- livsig[livsig$log2FoldChange < 0.00, ]
df_DEG <- read.csv("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/liv_df_DEG_v2.csv")
liv.res.df1.down_regulated <- df_DEG[which(df_DEG$log2FoldChange < 0), ]


## clusterProfiler
# gene ontology for liver BirA/Sox2Cre enriched genes over BirA 
library(clusterProfiler)
library(ggplot2)

livRNA.cc <- enrichGO(gene         = livsig_neg$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

livRNA.bp <- enrichGO(gene         = livsig_neg$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

livRNA.mf <- enrichGO(gene         = livsig_neg$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

livRNA.cc <- enrichGO(gene         = liv.res.df1.down_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

livRNA.bp <- enrichGO(gene         = liv.res.df1.down_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

livRNA.mf <- enrichGO(gene         = liv.res.df1.down_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)


# no enriched terms found 
livRNA.cc_plot <- dotplot(livRNA.cc, font.size = 16)
livRNA.bp_plot <- dotplot(livRNA.bp, font.size = 16)
livRNA.mf_plot <- dotplot(livRNA.mf, font.size = 16)

ggsave(filename = "GO_livRNA_neg_LG_cc_plot_v2.pdf", plot = livRNA.cc_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_livRNA_neg_LG_bp_plot_v2.pdf", plot = livRNA.bp_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_livRNA_neg_LG_mf_plot_v2.pdf", plot = livRNA.mf_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/liver/",  width = 16, height = 12, units = "in", dpi = "retina")

##################################################### DESeq Kidney BirA/Sox2Cre to BirA #######################################################
# make count table 
kid.counts <- data.frame("Gene" = kid92$Feature, "Kid92" = kid92$Count, "Kid93" = kid93$Count, "Kid94" = kid94$Count, "Kid95" = kid95$Count)

# make meta data table
kid.meta <- data.frame("id" = c("Kid92", "Kid93", "Kid94", "Kid95"), 
                       "genotype" = factor(c("BirA", "BirA", "BirA-Sox2Cre", "BirA-Sox2Cre"), levels = c("BirA-Sox2Cre", "BirA")), 
"tissue" = rep(c("kidney"), n = 4))


## deseq2 setup
kid.dds <- DESeqDataSetFromMatrix(countData = kid.counts, colData = kid.meta, design = ~genotype, tidy = TRUE)
saveRDS(kid.dds, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/kidney_BirA_over_BirA-Sox2Cre_v2.Rds", compress = FALSE)
kid.dds <- readRDS("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/kidney_BirA_over_BirA-Sox2Cre_v2.Rds")

# (Data has to be log transformed for that)
kid.dds_Rlog <- rlog(kid.dds, blind = FALSE)

plotPCA(kid.dds_Rlog, intgroup = "id", ntop = nrow(kid.dds_Rlog)) + theme_classic()
plotPCA(kid.dds_Rlog, intgroup = "genotype", ntop = nrow(kid.dds_Rlog)) + theme_classic()


# differential expression analysis
kid.dds <- DESeq(kid.dds)

# results table
kid.res <- results(kid.dds)
head(results(kid.dds, tidy=TRUE))

# summary 
summary(kid.res)


# sorted by p-value
kid.res <- kid.res[order(kid.res$padj), ]
head(kid.res)

# add in geneSymbols
kid.mapped <- kid.res
kid.mapped$SYMBOL <- rownames(kid.mapped)
# use library AnnotationDB
kid.mapped$symbol <- mapIds(org.Mm.eg.db, keys = row.names(kid.res), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
kidLog.mapped <- as.data.frame(kid.mapped@listData)

# filter based sig. p-value 
kidsig <- kidLog.mapped[kidLog.mapped$padj <= 0.05, ]
write.csv(kidsig, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/kidney_sig_BirA-BirA-Sox2Cre_v2.csv")

kidsig <- read.csv("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/kidney_sig_BirA-BirA-Sox2Cre_v2.csv")
kid.acc <- read.csv("/Users/amandameyer/Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/TR01_Sox2_LCMSMS/input/TR01_Kidney_mediannormalized_multimedian.csv")

##### from Louisa 
#NA omit
kid.res.df <- na.omit(as.data.frame(kid.res))
range(kid.res.df$baseMean)


mcols(kid.dds, use.names=TRUE)

#add normalized count of each sample
norm_counts_df <- as.data.frame(counts(kid.dds, normalized = T))
norm_counts_df <- norm_counts_df[which(row.names(norm_counts_df) %in% row.names(kid.res.df)), ]

COL_data <- as.data.frame(kid.dds@colData)
BirA <- COL_data[which(COL_data$tissue == "BirA"), 4]
BirA_Sox2Cre <- COL_data[which(COL_data$tissue == "BirA-Sox2Cre"), 4]


#Get mean for each condition
baseMean_BirA <- rowMeans(norm_counts_df[ , 1:2])
baseMean_BirA_Sox2Cre <- rowMeans(norm_counts_df[ , 3:4])


#combine df with normalized counts and mean for condition, sort df
norm_counts_df <- cbind(baseMean_BirA, baseMean_BirA_Sox2Cre, norm_counts_df)


kid.res.df1 <- cbind(kid.res.df, norm_counts_df)
head(kid.res.df1)

### Annotating Gene Symbol/Description with BiomaRt
kid.res.df1 <- kid.res.df1 %>% tibble::rownames_to_column("ensembl_gene_id") 

library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

kid_datasets.anno <- getBM(values = kid.res.df1$ensembl_gene_id, filters = "ensembl_gene_id", 
                           attributes = c("ensembl_gene_id", "entrezgene_id", "description", 'mgi_symbol'), mart = mart)

## use Louisa's nomenclature 
saveRDS(kid_datasets.anno, "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/kid_datasets_anno_v2.Rds", compress = FALSE)
anno <- readRDS("/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/kid_datasets_anno_v2.Rds")

kid.res.df1.anno <- merge(anno, kid.res.df1, by="ensembl_gene_id")
head(kid.res.df1.anno)

#Keep only genes with absolute logFC >= 2 and p value adj < 0.05
df_DEG <- kid.res.df1.anno %>% arrange(-log2FoldChange) %>% filter(padj < 0.05 & abs(log2FoldChange) >= 2)

#reorder df_DEG and save csv
#df_DEG <- df_DEG[ , c(1,4,3,9,10,6,7,8,2,5,11:length(colnames(df_DEG)))]

write.csv(df_DEG, file = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/kid_df_DEG_v2.csv", row.names = FALSE)

head(df_DEG, n = 20)

## Part 2-1: Volcano plot
## Create a column to indicate which genes to label
library(dplyr)
kid.res.df1.anno  <- mutate(kid.res.df1.anno , color = case_when(kid.res.df1.anno$log2FoldChange > 0 & kid.res.df1.anno$padj <0.05 ~ "Increased",
                                                                 kid.res.df1.anno$log2FoldChange < 0 & kid.res.df1.anno$padj <0.05 ~ "Decreased",
                                                                 kid.res.df1.anno$padj > 0.05 ~ "nonsignificant"))
kid.res.df1.anno <- kid.res.df1.anno %>% mutate(threshold = padj < 0.05) %>% arrange(padj) %>% mutate(volcanolabels = "")

#Top20
kid.res.df1.anno$volcanolabels[1:20] <- kid.res.df1.anno$mgi_symbol[1:20]

library(ggplot2)
library(ggrepel)
vol_p <-  ggplot(kid.res.df1.anno, aes(x = log2FoldChange, y = -log10(padj), color=color)) +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  geom_text_repel(aes(label = volcanolabels)) + scale_color_manual(name = "Directionality",
                                                                   values = c(Increased = "#D53E4F", Decreased = "#9ECAE1", nonsignificant = "grey50"))+
  xlab("log2(Fold Change)") + 
  ylab("-log10(Adj. p-value)") + theme_classic()
vol_p
ggsave(filename = "kid_DEG_volcano_plot_v2.pdf", plot = vol_p, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney", width = 6, height = 4, units = "in", dpi = "retina")

#"grey50", "#9ECAE1", "#C6DBEF", "#D53E4F"
## Part 2-2: Heatmap of top 15 up and top 15 downregulated genes
heatmap_genes <- df_DEG %>% arrange(padj) %>% top_n(n = -30, wt = padj)

# Gather 30 significant genes and make matrix
mat <- assay(kid.dds_Rlog)[heatmap_genes$ensembl_gene_id, ]
row.names(mat) <- heatmap_genes$mgi_symbol

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(kid.dds_Rlog)$Genotype), 
  row.names = colData(kid.dds_Rlog)$rownames
)

kid.dds_Rlog@colData$genotype

# Specify colors you want to annotate the columns by. Can add more groups like replicates here
ann_colors = list(Group = c(BirA = "grey50", BirA-Sox2Cre = "#9ECAE1"))

# Make Heatmap with pheatmap function.
## See more in documentation for customization
library(pheatmap)
plot_heatmap <- pheatmap(mat = mat, 
                         color = colorRampPalette(brewer.pal(9, "Blues"))(255), 
                         scale = "row", # Scale genes to Z-score (how many standard deviations)
                         annotation_colors = c("grey50", "#9ECAE1"),# Change the default colors of the annotations
                         fontsize = 6.5, # Make fonts smaller
                         cellwidth = 15, # Make the cells wider
                         show_colnames = T)



##Single gene plot
# Get gene with highest expression
top_gene <- kid.res.df1.anno$ensembl_gene_id[which.max(kid.res.df1.anno$log2FoldChange)]
plot_title <- kid.res.df1.anno$mgi_symbol[which(kid.res.df1.anno$ensembl_gene_id == top_gene)]

# Plot single gene
plotCounts(dds = kid.dds, 
           gene = top_gene, 
           intgroup = "genotype", 
           normalized = T, 
           transform = T, main =plot_title)



################################################################ GO analysis Kidney ###########################################################################
# BirA enriched genes over BirA/Sox2Cre
kidsig_pos <- kidsig[kidsig$log2FoldChange > 0.00, ]
kid.res.df1.up_regulated <- df_DEG[which(df_DEG$log2FoldChange > 0), ]


## clusterProfiler
# gene ontology for kidney BirA enriched genes over BirA/Sox2Cre 
library(clusterProfiler)
library(ggplot2)

kidRNA.cc <- enrichGO(gene         = kidsig_pos$SYMBOL,
                         OrgDb         = "org.Mm.eg.db",
                         keyType       = 'ENSEMBL',
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

kidRNA.bp <- enrichGO(gene         = kidsig_pos$SYMBOL,
                         OrgDb         = "org.Mm.eg.db",
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

kidRNA.mf <- enrichGO(gene         = kidsig_pos$SYMBOL,
                         OrgDb         = "org.Mm.eg.db",
                         keyType       = 'ENSEMBL',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

kidRNA.cc <- enrichGO(gene         = kid.res.df1.up_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

kidRNA.bp <- enrichGO(gene         = kid.res.df1.up_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

kidRNA.mf <- enrichGO(gene         = kid.res.df1.up_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)


# no enriched terms found 
kidRNA.cc_plot <- dotplot(kidRNA.cc, font.size = 16)
kidRNA.bp_plot <- dotplot(kidRNA.bp, font.size = 16)
kidRNA.mf_plot <- dotplot(kidRNA.mf, font.size = 16)

ggsave(filename = "GO_kidRNA_pos_LG_cc_plot_v2.pdf", plot = kidRNA.cc_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_kidRNA_pos_LG_bp_plot_v2.pdf", plot = kidRNA.bp_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_kidRNA_pos_LG_mf_plot_v2.pdf", plot = kidRNA.mf_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/",  width = 16, height = 12, units = "in", dpi = "retina")

##### negative 
# BirA/Sox2Cre enriched genes over BirA
kidsig_neg <- kidsig[kidsig$log2FoldChange > 0.00, ]
kid.res.df1.down_regulated <- df_DEG[which(df_DEG$log2FoldChange < 0), ]


## clusterProfiler
# gene ontology for kidney BirA/Sox2Cre enriched genes over BirA 
library(clusterProfiler)
library(ggplot2)

kidRNA.cc <- enrichGO(gene         = kidsig_neg$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

kidRNA.bp <- enrichGO(gene         = kidsig_neg$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

kidRNA.mf <- enrichGO(gene         = kidsig_neg$SYMBOL,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

kidRNA.cc <- enrichGO(gene         = kid.res.df1.down_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

kidRNA.bp <- enrichGO(gene         = kid.res.df1.down_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

kidRNA.mf <- enrichGO(gene         = kid.res.df1.down_regulated$ensembl_gene_id,
                      OrgDb         = "org.Mm.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)


# no enriched terms found 
kidRNA.cc_plot <- dotplot(kidRNA.cc, font.size = 16)
kidRNA.bp_plot <- dotplot(kidRNA.bp, font.size = 16)
kidRNA.mf_plot <- dotplot(kidRNA.mf, font.size = 16)

ggsave(filename = "GO_kidRNA_neg_LG_cc_plot_v2.pdf", plot = kidRNA.cc_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_kidRNA_neg_LG_bp_plot_v2.pdf", plot = kidRNA.bp_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_kidRNA_neg_LG_mf_plot_v2.pdf", plot = kidRNA.mf_plot, path = "/Volumes/LaCie 1/RNAseq/BirA_Sox2Cre/WashU_analysis/output/kidney/",  width = 16, height = 12, units = "in", dpi = "retina")

########################################################## DESeq by tissue type #############################################################
## deseq2 setup
allTT.dds <- DESeqDataSetFromMatrix(countData = allT.counts, colData = allT.meta, design = ~tissue, tidy = TRUE)

ddsLRT <- DESeq(dds1, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)

# brain, kidney, liver (order)
levels(allTT.dds$tissue)


# differential expression analysis 3 way 
allTT.dds <- DESeq(allTT.dds, test = "LRT", reduced = ~1)
resLRT <- results(allTT.dds)
resLRT.brk <- results(allTT.dds, contrast = c("tissue", "brain", "kidney"))
summary(resLRT)

results(allTT.dds, name = "tissue_brain_vs_kidney" , test="Wald")



# results table 
allTT.res <- results(allTT.dds)
head(results(allTT.dds, tidy=TRUE))

# summary 
summary(allTT.res)

# sorted by p-value
allTT.res <- allT.res[order(allTT.res$padj), ]
head(allTT.res)






library(biomaRt)
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

allT.values <- as.vector(allT.counts$Gene)

# query biomart
MM.results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"),
                 filters = "ensembl_transcript_id", values = allT.values,
                 mart = mart)

MM.results
#   ensembl_gene_id ensembl_transcript_id ensembl_peptide_id
# 1 ENSG00000163734       ENST00000296026    ENSP00000296026



mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")


library("AnnotationDbi")
library("org.Mm.eg.db")


allTT.mapped <- allTT.res

allTT.mapped$SYMBOL <- rownames(allTT.mapped)

allTT.mapped$symbol <- mapIds(org.Mm.eg.db, keys = row.names(allTT.res), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

log.mapped <- as.data.frame(allTT.mapped@listData)




