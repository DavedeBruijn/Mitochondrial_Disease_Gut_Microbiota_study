#Loading packages
library(phyloseq)
library(vegan)
library(dplyr)
library(rbiom)

# setting library paths
MIDLOCMicrobiome_abundance_biompath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/04_Cleaned/04_cleaned_MIDLOC.biom"

# import biom file into R
MIDLOCMicrobiome_abundance_biomdata <- as_rbiom(MIDLOCMicrobiome_abundance_biompath)
MIDLOCMicrobiome_abundance_biomdata

#renaming the ranks
colnames(MIDLOCMicrobiome_abundance_biomdata$taxonomy) <- c(".otu", "pathway", "genus", "species")
MIDLOCMicrobiome_abundance_biomdata

#normalize by the root since the colsum is not comparable
library(slam)
feature_table <- MIDLOCMicrobiome_abundance_biomdata$counts
col_sums <- col_sums(feature_table)
features_table_norm <- feature_table
features_table_norm$v <- features_table_norm$v/col_sums[feature_table$j]
col_sums(features_table_norm)
MIDLOCMicrobiome_abundance_biomdata$counts <- features_table_norm

#check the otu_table 
otu_table <- as.matrix(MIDLOCMicrobiome_abundance_biomdata$counts)
colSums(otu_table)
View(otu_table)


#alpha diversity (Shannon)
ad_MIDLOCMicrobiome_abundance_shannon  <- adiv_table(MIDLOCMicrobiome_abundance_biomdata, adiv = "Shannon")
View(ad_MIDLOCMicrobiome_abundance_shannon)

#alpha diversity (Simpson)
ad_MIDLOCMicrobiome_abundance_simpson  <- adiv_table(MIDLOCMicrobiome_abundance_biomdata, adiv = "Simpson")
View(ad_MIDLOCMicrobiome_abundance_simpson)

#alpha diversity (All analysis)
ad_MIDLOCMicrobiome_abundance <- adiv_table(MIDLOCMicrobiome_abundance_biomdata, adiv = ".all")
View(ad_MIDLOCMicrobiome_abundance)

#boxplot alpha diversity and statistical test (Can change the stat.by)(Unpaired analysis)
adiv_boxplot_microbiome_abundance <- adiv_boxplot(MIDLOCMicrobiome_abundance_biomdata, layers = "vdp", adiv = "Shannon", stat.by = "group", x = 'group')
adiv_boxplot_microbiome_abundance

#the statistical test for the alpha diversity (However we use paired wilcoxon test so we can't use this)
stats <- adiv_stats(MIDLOCMicrobiome_abundance_biomdata, stat.by = "group", adiv = "Shannon", test = "Wilcox")
stats

#creating a nicer plot
library(ggplot2)
adiv_microbial_abundance <- ggplot(ad_MIDLOCMicrobiome_abundance_shannon, aes(x = group, y = .diversity, color = group))+
  geom_violin(fill = NA, linewidth = 1.3)+
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9)+
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.75)+
  scale_color_manual(values = c("Control" = "#00BFC4", "MD" = "#F8766D", "T1D" = "#BA72FF"))+
  theme_minimal()+
  labs(y = "Shannon Diversity (Metacyc Pathways)")+
  annotate("segment", x = 1, xend = 3, y = 3, yend = 3, linewidth = 0.5, color = "black") +
  annotate("segment", x = 1, xend = 1, y = 2.95, yend = 3.05, linewidth = 0.5, color = "black") +
  annotate("segment", x = 3, xend = 3, y = 2.95, yend = 3.05, linewidth = 0.5, color = "black") +
  annotate("text", x = 2, y = 3.05, label = "p = 0.004*", size = 4)
adiv_microbial_abundance

#Saving the shannon diversity
library(readr)
write_tsv(ad_MIDLOCMicrobiome_abundance_shannon, "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/plots/adiv_functional.tsv")