#Loading packages
library(phyloseq)
library(vegan)
library(dplyr)
library(rbiom)

# setting library paths
MIDLOCMicrobiome_biompath <- 

# import biom file into R
MIDLOCMicrobiome_biomdata <- as_rbiom(MIDLOCMicrobiome_biompath)
MIDLOCMicrobiome_biomdata

#changing type to group
names(MIDLOCMicrobiome_biomdata$metadata)[2] <- "group"

#normalize by the root since certain count were removed
features_table <- MIDLOCMicrobiome_biomdata$counts
col_sums <- col_sums(features_table)
features_table_norm <- features_table
features_table_norm$v <- features_table_norm$v/col_sums[features_table$j]
col_sums(features_table_norm)
MIDLOCMicrobiome_biomdata$counts <- features_table_norm

#making otu_table
otu_table <- as.matrix(MIDLOCMicrobiome_biomdata$counts)

#alpha diversity (Shannon)
ad_MIDLOCMicrobiome_shannon  <- adiv_table(MIDLOCMicrobiome_biomdata, adiv = "Shannon")
View(ad_MIDLOCMicrobiome_shannon)

#alpha diversity (Simpson)
ad_MIDLOCMicrobiome_simpson  <- adiv_table(MIDLOCMicrobiome_biomdata, adiv = "Simpson")
View(ad_MIDLOCMicrobiome_simpson)

#alpha diversity (All analysis)
ad_MIDLOCMicrobiome <- adiv_table(MIDLOCMicrobiome_biomdata, adiv = ".all")
View(ad_MIDLOCMicrobiome)

#boxplot alpha diversity and statistical test (Can change the stat.by)(Unpaired analysis)
adiv_boxplot_microbiome <- adiv_boxplot(MIDLOCMicrobiome_biomdata, layers = "vdp", adiv = "Shannon", stat.by = "type", x = 'type')
adiv_boxplot_microbiome

#the statistical test for the alpha diversity (However we use paired wilcoxon test so we can't use this)
adiv_stats(MIDLOCMicrobiome_biomdata, stat.by = "type", adiv = "Shannon", test = "Wilcox")

#creating a nicer plot
library(ggplot2)
adiv_microbial <- ggplot(ad_MIDLOCMicrobiome_shannon, aes(x = type, y = .diversity, color = type))+
  geom_violin(fill = NA, linewidth = 1.3)+
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9)+
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.75)+
  theme_minimal()+
  labs(y = "Shannon Diversity (Microbial)")+
  ggtitle("Alpha Diversity of patients with MD compared to healthy controls")+
  annotate("segment", x = 1, xend = 2, y = 5.25, yend = 5.25, linewidth = 0.5, color = "black") +
  annotate("segment", x = 1, xend = 1, y = 5.15, yend = 5.35, linewidth = 0.5, color = "black") +
  annotate("segment", x = 2, xend = 2, y = 5.15, yend = 5.35, linewidth = 0.5, color = "black") +
  annotate("text", x = 1.5, y = 5.35, label = "p = 0.0603", size = 4)
adiv_microbial

#creating plot for the article
adiv_microbial_art <- ggplot(ad_MIDLOCMicrobiome_shannon, aes(x = group, y = .diversity, color = group))+
  geom_violin(fill = NA, linewidth = 1.3)+
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9)+
  scale_color_manual(values = c("Control" = "#00BFC4", "MD" = "#F8766D", "T1D" = "#BA72FF"))+
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.75)+
  theme_minimal()+
  labs(y = "Shannon Diversity (Microbial)")
adiv_microbial_art

#adding the functional Alpha diversity to combine it into one graph
#adding functional alpha diversity to this project so I could combine it
library(readr)
ad_MIDLOCMicrobiome_shannon_functional <- read_tsv("/adiv_functional.tsv")

adiv_microbial_art_func <- ggplot(ad_MIDLOCMicrobiome_shannon_functional, aes(x = group, y = .diversity, color = group))+
  geom_violin(fill = NA, linewidth = 1.3)+
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9)+
  scale_color_manual(values = c("Control" = "#00BFC4", "MD" = "#F8766D", "T1D" = "#BA72FF"))+
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.75)+
  theme_minimal()+
  labs(y = "Shannon Diversity (Metacyc Pathways)")+
  annotate("segment", x = 1, xend = 3, y = 2.95, yend = 2.95, linewidth = 0.5, color = "black") +
  annotate("segment", x = 1, xend = 1, y = 2.9, yend = 3, linewidth = 0.5, color = "black") +
  annotate("segment", x = 3, xend = 3, y = 2.9, yend = 3, linewidth = 0.5, color = "black") +
  annotate("text", x = 2, y = 3.05, label = "p = 0.004*", size = 4)
adiv_microbial_art_func

#combining the plot
library(patchwork)
adiv_complete <- adiv_microbial_art / adiv_microbial_art_func
adiv_complete + plot_annotation(
  tag_levels = "A",
  theme = theme(plot.title = element_text(hjust = 0.5))
)
