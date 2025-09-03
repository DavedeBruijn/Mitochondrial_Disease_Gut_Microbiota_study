#Loading packages
library(phyloseq)
library(vegan)
library(dplyr)
library(rbiom)

# setting library paths
MIDLOCMicrobiome_biompath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/04_cleaned/04_cleaned_MIDLOC.biom"

# import biom file into R
MIDLOCMicrobiome_biomdata <- as_rbiom(MIDLOCMicrobiome_biompath)
MIDLOCMicrobiome_biomdata

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

