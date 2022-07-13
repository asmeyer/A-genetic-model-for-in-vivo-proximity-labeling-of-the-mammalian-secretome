# 2021-03-22
# TR01: BirA/Sox2 Serum 
# ASM 

# last updated: 2022-01-20

# usage: script to make volcano plot of enriched proteins in a 6plex TMT mass spec experiment. 
      # Note:" ASM = serum plex1 and Ilia = serum plex2


####################################################################### begin ##########################################################################
library(dplyr)
library(tidyr)
library(tidyverse)

# data from two separate runs: sample prep by Ilia (I) or Amanda (A)


# will start with ASM data to make pipeline, then will run same analysis for Ilia's data 
# ASM not normalized multi-median group comparison for log2 FC data
serum_Am <- read.csv("TR01_Serum_AmandaPlex_Multimedian_notNormalized_Unfiltered_correctOrder.csv")
head(serum_Am)
serum_Am 

serum_Am.filt <- serum_Am[serum_Am$Amanda_newDatabase.unique_peptides >= 2, ]
str(serum_Am.filt)
str(serum_Am.filt$Amanda_newDatabase.pct_cov)
length(serum_Am.filt$Amanda_newDatabase.unique_peptides)
length(serum_Am$Amanda_newDatabase.unique_peptides)

# ASM NORMALIZED multi-median group comparison for log2 FC data (will not use, normalization causes issue; will show)
serum_Am.norm <- read.csv("TR01_Serum_AmandaPlex_Multimedian_Normalized_Unfiltered_correctOrder.csv") # serum_A 
head(serum_Am.norm)
serum_Am.norm 

# Ilia not normalized multi-median group comparison for log2 FC data
serum_Im <- read.csv("TR01_Serum_IliaPlex_Multimedian_notNormalized_Unfiltered_correctOrder.csv") # serum_A 
head(serum_Im)
serum_Im
# Ilia NORMALIZED multi-median group comparison for log2 FC data (will not use, normalization causes issue; will show)
serum_Im.norm <- read.csv("TR01_Serum_IliaPlex_Multimedian_Normalized_Unfiltered_correctOrder.csv") # serum_A 
head(serum_Im.norm)
serum_Im.norm



####################################################################### label points ##########################################################################
# label samples based on enrichment identity 
# ASM notNorm
serum_Am$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
serum_Am$diffexpressed[serum_Am$logFC.xBirA.Cre.over.NoBirA > 1.0 & serum_Am$P.Value.xBirA.Cre.over.NoBirA < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
serum_Am$diffexpressed[serum_Am$logFC.xBirA.Cre.over.NoBirA < -1.0 & serum_Am$P.Value.xBirA.Cre.over.NoBirA < 0.05] <- "DOWN"
serum_Am$diffexpressed[serum_Am$logFC.xBirA.Cre.over.NoBirA > 1.0 & serum_Am$P.Value.xBirA.Cre.over.NoBirA > 0.05] <- "NS_UP"
sum(serum_Am$diffexpressed == "UP")

# label samples based on enrichment identity 
# ASM notNorm filtered 
serum_Am.filt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
serum_Am.filt$diffexpressed[serum_Am.filt$logFC.xBirA.Cre.over.NoBirA > 1.0 & serum_Am.filt$P.Value.xBirA.Cre.over.NoBirA < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
serum_Am.filt$diffexpressed[serum_Am.filt$logFC.xBirA.Cre.over.NoBirA < -1.0 & serum_Am.filt$P.Value.xBirA.Cre.over.NoBirA < 0.05] <- "DOWN"
serum_Am.filt$diffexpressed[serum_Am.filt$logFC.xBirA.Cre.over.NoBirA > 1.0 & serum_Am.filt$P.Value.xBirA.Cre.over.NoBirA > 0.05] <- "NS_UP"
sum(serum_Am.filt$diffexpressed == "UP")

# ASM Norm
serum_Am.norm$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
serum_Am.norm$diffexpressed[serum_Am.norm$logFC.xBirA.Cre.over.NoBirA > 1.0 & serum_Am.norm$P.Value.xBirA.Cre.over.NoBirA < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
serum_Am.norm$diffexpressed[serum_Am.norm$logFC.xBirA.Cre.over.NoBirA < -1.0 & serum_Am.norm$P.Value.xBirA.Cre.over.NoBirA < 0.05] <- "DOWN"
serum_Am.norm$diffexpressed[serum_Am.norm$logFC.xBirA.Cre.over.NoBirA > 1.0 & serum_Am.norm$P.Value.xBirA.Cre.over.NoBirA > 0.05] <- "NS_UP"


# Ilia notNorm
serum_Im$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
serum_Im$diffexpressed[serum_Im$logFC.xBirA.Cre.over.NoBirA > 1.0 & serum_Im$P.Value.xBirA.Cre.over.NoBirA < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
serum_Im$diffexpressed[serum_Im$logFC.xBirA.Cre.over.NoBirA < -1.0 & serum_Im$P.Value.xBirA.Cre.over.NoBirA < 0.05] <- "DOWN"
serum_Im$diffexpressed[serum_Im$logFC.xBirA.Cre.over.NoBirA > 1.0 & serum_Im$P.Value.xBirA.Cre.over.NoBirA > 0.05] <- "NS_UP"
sum(serum_Im$diffexpressed == "UP")

# Ilia Norm
serum_Im.norm$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
serum_Im.norm$diffexpressed[serum_Im.norm$logFC.xBirA.Cre.over.NoBirA > 1.0 & serum_Im.norm$P.Value.xBirA.Cre.over.NoBirA < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
serum_Im.norm$diffexpressed[serum_Im.norm$logFC.xBirA.Cre.over.NoBirA < -1.0 & serum_Im.norm$P.Value.xBirA.Cre.over.NoBirA < 0.05] <- "DOWN"
serum_Im.norm$diffexpressed[serum_Im.norm$logFC.xBirA.Cre.over.NoBirA > 1.0 & serum_Im.norm$P.Value.xBirA.Cre.over.NoBirA > 0.05] <- "NS_UP"

################################################################## ASM Volcano #########################################################################
# actual plot to use
# plots plots plots 
plot1 <- ggplot(serum_Im, aes(x = logFC.xBirA.Cre.over.NoBirA, y = Log.P.Value.xBirA.Cre.over.NoBirA, col = diffexpressed, label = geneSymbol)) + 
  geom_point(alpha = 0.60) +
  #scale_color_manual(labels = c("NS", "Log2 FC", "p-value & log2 FC"), values = c("grey50", "#9ECAE1", "#D53E4F") +
  theme_classic(base_size = 30) + 
  #geom_text() +  
  geom_vline(xintercept = c(1.0), col = "black") +
  geom_hline(yintercept = -10*log10(0.05), col = "black")

plot1 <- plot1 + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(y = "-10log"[10]~"(p-value)") +
  labs(x =  "log"[2]~"fold change") + 
  labs(color = "") + 
  scale_color_manual(labels = c("log"[2]~"FC", "NS", "log"[2]~"FC", "p-value & log"[2]~"FC"), values = c("#6BAED6", "grey50", "#9ECAE1", "#D53E4F"))
plot1



plot1 <- ggplot(serum_Am.filt, aes(x = logFC.xBirA.Cre.over.NoBirA, y = Log.P.Value.xBirA.Cre.over.NoBirA, col = diffexpressed, label = geneSymbol)) + 
  geom_point(alpha = 0.60) +
  #scale_color_manual(labels = c("NS", "Log2 FC", "p-value & log2 FC"), values = c("grey50", "#9ECAE1", "#D53E4F") +
  theme_classic(base_size = 30) + 
  #geom_text() +  
  geom_vline(xintercept = c(1.0), col = "black") +
  geom_hline(yintercept = -10*log10(0.05), col = "black")

plot1 <- plot1 + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(y = "-10log"[10]~"(p-value)") +
  labs(x =  "log"[2]~"fold change") + 
  labs(color = "") + 
  scale_color_manual(labels = c("p-value & -log"[2]~"FC", "NS", "log"[2]~"FC", "p-value & log"[2]~"FC"), values = c("#6BAED6", "grey50", "#9ECAE1", "#D53E4F"))
plot1




# ASM notNorm
ggsave(filename = "BirA-Sox2Cre_Serum_LCMSMS_multimed_notNorm_ASM_volcano.pdf", plot = plot1, path = "./", 
       width = 14, height = 12, units = "in", dpi = "retina")

# filt
ggsave(filename = "BirA-Sox2Cre_Serum_LCMSMS_multimed_notNorm-filt_UniquePeps_ASM_volcano.pdf", plot = plot1, path = "./", 
       width = 14, height = 12, units = "in", dpi = "retina")

# ASM Norm
ggsave(filename = "BirA-Sox2Cre_Serum_LCMSMS_multimed_Norm_ASM_volcano.pdf", plot = plot1, path = "./", 
       width = 14, height = 12, units = "in", dpi = "retina")


# Ilia notNorm
ggsave(filename = "BirA-Sox2Cre_Serum_LCMSMS_multimed_notNorm_Ilia_volcano.pdf", plot = plot1, path = "./", 
       width = 14, height = 12, units = "in", dpi = "retina")

# Ilia Norm
ggsave(filename = "BirA-Sox2Cre_Serum_LCMSMS_multimed_Norm_Ilia_volcano.pdf", plot = plot1, path = "./", 
       width = 14, height = 12, units = "in", dpi = "retina")
