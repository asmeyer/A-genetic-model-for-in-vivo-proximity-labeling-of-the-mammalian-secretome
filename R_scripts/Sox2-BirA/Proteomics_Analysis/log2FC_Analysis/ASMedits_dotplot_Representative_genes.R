###################Dotchart to show representative genes###################
# import data
library("readxl")
rep_genes <- read_excel("Representative_unique_genes.xlsx")
rep_genes
# apply ggdotchart
library(ggpubr)
ggdotchart(rep_genes, x="genes",y="BirA;Sox2-Cre/BirA fold change",
           groups = "Tissue_type", color = "Tissue_type",
           palette = c('#999999','#E69F00','#56B4E9','#D4596F'),
           rotate = F,position = "dodge",
           ggtheme = theme_bw(),
           y.text.col = TRUE )

# make plot 
plot1 <- ggplot(rep_genes, aes(genes, `BirA;Sox2-Cre/BirA fold change`, color = Tissue_type)) + 
  geom_point(size = 5, alpha = 0.75) +
  facet_wrap(~ Tissue_type, scales = "free_x", ncol = 4) + 
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 45, hjust = 0.8, vjust = 0.75, face = "plain")) +
  scale_color_manual(values = c("#257DCF", "#42853D", "#AA0A3C", "#767679")) +
  ggtitle("Log"[2]~"fold change BirA-Sox2Cre/BirA") + theme(plot.title = element_text(hjust = 0.5)) +
  labs(caption = "") +
  labs(y = "log"[2]~"fold change") +
  labs(x = "Protein") + 
  labs(fill = "Tissue")

# print plot 
plot1

# save plot 
ggsave(filename = "Point_plot_log2_proteins_kid-liv-brain-ser.pdf", plot = plot1, path = "Figures/", width = 14, height = 12, units = "in", dpi = "retina")



# make plot 
plot1 <- ggplot(rep_genes, aes(genes, `BirA;Sox2-Cre/BirA fold change`, fill = genes)) + 
  geom_bar(stat = "identity", position=position_dodge(), alpha = 0.65) +
  facet_wrap(~ Tissue_type, scales = "free_x", ncol = 4) + 
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 45, hjust = 0.8, vjust = 0.75, face = "plain")) +
  scale_fill_manual(values = c("#257DCF", "#42853D", "#AA0A3C", "#767679")) +
  ggtitle("Log"[2]~"fold change BirA-Sox2Cre/BirA") + theme(plot.title = element_text(hjust = 0.5)) +
  labs(caption = "") +
  labs(y = "log"[2]~"fold change") +
  labs(x = "Protein") + 
  labs(fill = "Tissue")

# print plot 
plot1

# save plot 
ggsave(filename = "Bar_plot_log2_proteins_kid-liv-brain-ser.pdf", plot = plot1, path = "Figures/", width = 14, height = 12, units = "in", dpi = "retina")


####### by protein 
# make plot 
# make plot 
plot1 <- ggplot(rep_genes, aes(Tissue_type, `BirA;Sox2-Cre/BirA fold change`, color = Tissue_type)) + 
  geom_point(size = 5, alpha = 0.75) +
  facet_wrap(~ genes, scales = "free_x", ncol = 4) + 
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 45, hjust = 0.8, vjust = 0.75, face = "plain")) +
  scale_color_manual(values = c("#257DCF", "#42853D", "#AA0A3C", "#767679")) +
  ggtitle("Log"[2]~"fold change BirA-Sox2Cre/BirA") + theme(plot.title = element_text(hjust = 0.5)) +
  labs(caption = "") +
  labs(y = "log"[2]~"fold change") +
  labs(x = "Protein") + 
  labs(fill = "Tissue")

# print plot 
plot1

# save plot 
ggsave(filename = "Point_plot_log2_proteins_kid-liv-brain-ser_by_Tissue.pdf", plot = plot1, path = "Figures/", width = 14, height = 12, units = "in", dpi = "retina")


# make plot 
plot1 <- ggplot(rep_genes, aes(genes, `BirA;Sox2-Cre/BirA fold change`, fill = genes, color = genes)) + 
  geom_bar(stat = "identity", position=position_dodge(), alpha = 0.65) +
  facet_wrap(~ Tissue_type, scales = "free_x", ncol = 4) + 
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 45, hjust = 0.8, vjust = 0.75, face = "plain")) +
  scale_color_manual(values = c("black", "black", "black", "black")) +
  scale_fill_manual(values = c("white", "grey80", "grey50", "#11161E")) +
  ggtitle("Log"[2]~"fold change BirA-Sox2Cre/BirA") + theme(plot.title = element_text(hjust = 0.5)) +
  labs(caption = "") +
  labs(y = "log"[2]~"fold change") +
  labs(x = "Protein") + 
  labs(fill = "Tissue")

# print plot 
plot1

# save plot 
ggsave(filename = "Bar_plot_log2_proteins_kid-liv-brain-ser_greys.pdf", plot = plot1, path = "Figures/", width = 14, height = 12, units = "in", dpi = "retina")




# make plot 
plot1 <- ggplot(rep_genes, aes(genes, `BirA;Sox2-Cre/BirA fold change`, fill = genes)) + 
  geom_bar(stat = "identity", position=position_dodge(), alpha = 0.65) +
  facet_wrap(~ Tissue_type, scales = "free_x", ncol = 4) + 
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 45, hjust = 0.8, vjust = 0.75, face = "plain")) +
  scale_fill_manual(values = c("#257DCF", "#42853D", "#AA0A3C", "#767679")) +
  ggtitle("Log"[2]~"fold change BirA-Sox2Cre/BirA") + theme(plot.title = element_text(hjust = 0.5)) +
  labs(caption = "") +
  labs(y = "log"[2]~"fold change") +
  labs(x = "Protein") + 
  labs(fill = "Tissue")

# print plot 
plot1

# save plot 
ggsave(filename = "Bar_plot_log2_proteins_kid-liv-brain-ser.pdf", plot = plot1, path = "Figures/", width = 14, height = 12, units = "in", dpi = "retina")

