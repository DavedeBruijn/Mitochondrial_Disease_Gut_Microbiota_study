#Loading packages
library(rbiom)
library(phyloseq)
library(vegan)
library(dplyr)
library(cowplot)
library(OmicFlow)
library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(slam)

# setting library paths
MIDLOCMicrobiome_biompath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/04_cleaned/04_cleaned_MIDLOC.biom"

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

#create an otu_table 
otu_table <- as.matrix(MIDLOCMicrobiome_biomdata$counts)

#Differential abundance analysis
Diff_Family_table <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = -3, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
Diff_genus_table <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = -2, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
Diff_species_table <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = -1, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_Family_table)
View(Diff_genus_table)
View(Diff_species_table)

#create a data.table 
genus_table <- taxa_table(MIDLOCMicrobiome_biomdata, rank = -2, taxa = 0.00000001, unc = "drop")
View(genus_table)
genus_list <- unique(genus_table$.taxa)
View(genus_list)

genus_table <- as.data.table(genus_table)
genus_datatable <- dcast(genus_table, .taxa ~ .sample, value.var = ".abundance")
View(genus_datatable)

#Adding the label for the conditions
setnames(genus_datatable, 
         old = names(genus_datatable)[2:31], 
         new = paste0(names(genus_datatable)[2:31], "_MD"))

# Rename columns 32 to 91 with _C suffix
setnames(genus_datatable, 
         old = names(genus_datatable)[32:91], 
         new = paste0(names(genus_datatable)[32:91], "_C"))
# Rename columns 92 to 151 with _C suffix
setnames(genus_datatable, 
         old = names(genus_datatable)[92:151], 
         new = paste0(names(genus_datatable)[92:151], "_T1D"))

#using the foldchange
DAA_genus <- foldchange(data = genus_datatable, feature_rank = ".taxa", condition_A = c("_MD", "_MD"), condition_B = c("_C", "_T1D") ,paired = FALSE, condition_labels = names(genus_datatable)[-1])
View(DAA_genus)

DAA_genus$FDR_1 <-p.adjust(DAA_genus$pvalue_1 , method = "fdr")
DAA_genus$FDR_2 <-p.adjust(DAA_genus$pvalue_2 , method = "fdr")

#checking with OmicFlow
#changing type to group
metadata <- as.data.table(MIDLOCMicrobiome_biomdata$metadata)
countData <- as.matrix(MIDLOCMicrobiome_biomdata$counts)
featureData <- as.data.frame(MIDLOCMicrobiome_biomdata$taxonomy)

colnames(metadata)[1]<- "SAMPLE_ID"
colnames(featureData)[1] <- "FEATURE_ID"

names(metadata)[2] <- "CONTRAST_Group"
metadata$CONTRAST_Group <- as.character(metadata$CONTRAST_Group)

taxa <- metagenomics$new(
  metaData = metadata,
  countData = countData,
  featureData = featureData
)
res <- taxa$DFE(
  feature_rank = "Genus",
  paired = FALSE,
  normalize = FALSE,
  condition.group = "CONTRAST_Group",
  condition_A = "MD",
  condition_B = "Control"
)

View(res$data)

#calculating the relative abudance for the volcano plot
rel_abun <- as.data.frame(res$data)
colnames(rel_abun)[1] <- ".taxa"

DAA_genus <- merge(DAA_genus, rel_abun[,c(".taxa", "rel_abun")], by = ".taxa", all.x = TRUE)
DAA_genus$rel_abun <- (DAA_genus$rel_abun/sum(DAA_genus$rel_abun))*100

#making a volcano plot from the different abundance analysis Genus
volcanoplot_results <- DAA_genus %>%
  mutate(log10_adj.pval = -log10(FDR_1),
         significant = FDR_1 < 0.05)

#removing Microbacterium from the plot to add manually for visability
Microbacterium <- volcanoplot_results[1013]
Sphingomonas <- volcanoplot_results[1070]
volcanoplot_results <- volcanoplot_results[-c(1013,1070)]

top_hits <- volcanoplot_results %>%
  filter(FDR_1 < 0.05) %>%
  arrange(FDR_1) 
top_hits_noGGB <- top_hits[!grepl("^GGB", top_hits$.taxa), ]

hline_df <- data.frame(yintercept = -log10(0.05), Label = "FDR = 0.05")

DAA_genus_plot <- ggplot(volcanoplot_results, aes(x = Log2FC_1, y = log10_adj.pval, colour = Log2FC_1))+
  geom_point(alpha = 1, size = 1) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
  theme_minimal()+
  labs(x = "log2 Fold Change(MD/Control)",
       y = "-log10(FDR)",
       size = "Relative Abundance (%)",
       colour = "Fold Change",
       title = "Genus")+
  geom_hline(data = hline_df, aes(yintercept = yintercept), color = "blue", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = -21, y = -log10(0.035), label = "FDR < 0.05", color = "black", size = 3)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = top_hits, aes(color = Log2FC_1, size = rel_abun), alpha = 1)+
  scale_size_continuous(range = c(0.1, 4), trans = "log10",breaks = c(1, 0.1, 0.01, 0.001), labels = c("1", "0.1", "0.01", "0.001"))+
  geom_text_repel(data = top_hits_noGGB, aes(x = Log2FC_1, y = log10_adj.pval, label = .taxa), color = "black", size = 2.5, force = 3.5, max.overlaps = 20)
DAA_genus_plot

DAA_genus_plot + plot_annotation(
  title = "Difference in relative abundance of genus in patients with MD and healthy controls",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)
legend <- get_legend(DAA_genus_plot + theme(legend.position="bottom", legend.direction = "horizontal"))
DAA_genus_plot <- DAA_genus_plot + theme(legend.position="none")

#making a species foldchange volcanoplot
#using OmicFlow to get the foldchange for the volcanoplot
species_table <- taxa_table(MIDLOCMicrobiome_biomdata, rank = -1, taxa = 0.00000001, unc = "drop")
View(species_table)
species_list <- unique(species_table$.taxa)
View(species_list)

species_table <- as.data.table(species_table)
species_datatable <- dcast(species_table, .taxa ~ .sample, value.var = ".abundance")
View(species_datatable)

#Adding the label for the conditions
setnames(species_datatable, 
         old = names(species_datatable)[2:31], 
         new = paste0(names(species_datatable)[2:31], "_MD"))

# Rename columns 32 to 91 with _C suffix
setnames(species_datatable, 
         old = names(species_datatable)[32:91], 
         new = paste0(names(species_datatable)[32:91], "_C"))

# Rename columns 92 to 151 with _T1D suffix
setnames(species_datatable, 
         old = names(species_datatable)[92:151], 
         new = paste0(names(species_datatable)[92:151], "_T1D"))

#using the foldchange
DAA_species <- foldchange(data = species_datatable, feature_rank = ".taxa", condition_A = c("_MD", "_MD"), condition_B = c("_C", "_T1D") ,paired = FALSE, condition_labels = names(species_datatable)[-1])
View(DAA_species)

DAA_species$FDR_1 <-p.adjust(DAA_species$pvalue_1 , method = "fdr")
DAA_species$FDR_2 <-p.adjust(DAA_species$pvalue_2 , method = "fdr")

#calculating the relative abudance for the volcano plot
res <- taxa$DFE(
  feature_rank = "Species",
  paired = FALSE,
  normalize = FALSE,
  condition.group = "CONTRAST_Group",
  condition_A = "MD",
  condition_B = "Control"
)

View(res$data)

#calculating the relative abudance for the volcano plot
rel_abun_species <- as.data.frame(res$data)
colnames(rel_abun_species)[1] <- ".taxa"

DAA_species <- merge(DAA_species, rel_abun_species[,c(".taxa", "rel_abun")], by = ".taxa", all.x = TRUE)
DAA_species$rel_abun <- (DAA_species$rel_abun/sum(DAA_species$rel_abun))*100


#making a volcano plot from the different abundance analysis Genus
volcanoplot_results_species <- DAA_species %>%
  mutate(log10_adj.pval = -log10(FDR_1),
         significant = FDR_1 < 0.05)

#Creating point outside the plot for visibility
Microbacterium_sp_T32 <- volcanoplot_results_species[1435,]
Sphingomonas_sp_FARSPH <- volcanoplot_results_species[1610,]
volcanoplot_results_species <- volcanoplot_results_species[-c(1435,1610),]

top_hits_species <- volcanoplot_results_species %>%
  filter(FDR_1 < 0.05) %>%
  arrange(FDR_1) 
top_hits_species_noGGB <- top_hits_species[!grepl("^GGB", top_hits_species$.taxa), ]

hline_df <- data.frame(yintercept = -log10(0.05), Label = "FDR = 0.05")

DAA_species_plot <- ggplot(volcanoplot_results_species, aes(x = Log2FC_1, y = log10_adj.pval, colour = Log2FC_1))+
  geom_point(alpha = 1, size = 1) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
  theme_minimal()+
  labs(x = "log2 Fold Change(MD/Control)",
       y = "-log10(FDR)",
       size = "Relative Abundance",
       colour = "Fold Change",
       title = "Species")+
  geom_hline(data = hline_df, aes(yintercept = yintercept), color = "blue", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = -21, y = -log10(0.035), label = "FDR < 0.05", color = "black", size = 3)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = top_hits_species, aes(color = Log2FC_1, size = rel_abun), alpha = 1)+
  scale_size_continuous(range = c(0.1, 4), trans = "log10",breaks = c(1, 0.1, 0.01, 0.001), labels = c("1", "0.1", "0.01", "0.001"))+
  geom_text_repel(data = top_hits_species_noGGB, aes(x = Log2FC_1, y = log10_adj.pval, label = .taxa), color = "black", size = 2.6, force = 3.5, max.overlaps = 16)

DAA_species_plot + plot_annotation(
  title = "Difference in relative abundance of species in patients with MD and healthy controls",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)
DAA_species_plot <- DAA_species_plot + theme(legend.position = "none")

#combining both volcano plots for the article
DAA_complete <- DAA_genus_plot/DAA_species_plot

DAA_complete <- DAA_complete + plot_annotation(
  tag_levels = "A",
)

DAA_complete <- cowplot::plot_grid(DAA_complete, legend, 
                                   nrow = 2, rel_heights = c(1, 0.1))
DAA_complete

#making the tables for the article
#making the Control - MD comparison
MIDLOC_C_MD <- subset(MIDLOCMicrobiome_biomdata, MIDLOCMicrobiome_biomdata$metadata$group != "T1D")
MIDLOC_C_MD

Diff_genus_table_C_MD <- taxa_stats(MIDLOC_C_MD, rank = -2, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_genus_table_C_MD)
Diff_species_table_C_MD <- taxa_stats(MIDLOC_C_MD, rank = -1, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_species_table_C_MD)

#Adding the foldchange to tables
Diff_genus_table_C_MD <- merge(Diff_genus_table_C_MD, DAA_genus, by = ".taxa")

#making the differential abundance tables for the article
genus_C_MD_article <- Diff_genus_table_C_MD[,-c(4,5,6,10, 13, 14,16)]
names(genus_C_MD_article)[7] <- "Log2FC"
names(genus_C_MD_article)[8] <- "Pvalue"
names(genus_C_MD_article)[9] <- "FDR"

#subsetting only the significant
genus_C_MD_article <- genus_C_MD_article %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)
View(genus_C_MD_article)

#Adding the foldchange to tables
Diff_species_table_C_MD <- merge(Diff_species_table_C_MD, DAA_species, by = ".taxa")

#making the differential abundance tables for the article
Species_C_MD_article <- Diff_species_table_C_MD[,-c(4,5,6,10, 13, 14,16)]
names(Species_C_MD_article)[7] <- "Log2FC"
names(Species_C_MD_article)[8] <- "Pvalue"
names(Species_C_MD_article)[9] <- "FDR"

#subsetting only the significant
Species_C_MD_article <- Species_C_MD_article %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)
View(Species_C_MD_article)

#creating the PCoA plot with only MD and T1D
MIDLOC_MD_T1D <- subset(MIDLOCMicrobiome_biomdata, MIDLOCMicrobiome_biomdata$metadata$group != "Control")
MIDLOC_MD_T1D

#swap the levels of the groups so that it becomes: MD control
MIDLOC_MD_T1D$metadata$group <- factor(MIDLOC_MD_T1D$metadata$group, levels = rev(levels(MIDLOC_MD_T1D$metadata$group)))


#Differential abundance analysis table MD_T1D
Diff_genus_table_MD_T1D <- taxa_stats(MIDLOC_MD_T1D, rank = -2, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
Diff_species_table_MD_T1D <- taxa_stats(MIDLOC_MD_T1D, rank = -1, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_genus_table_MD_T1D)
View(Diff_species_table_MD_T1D)

#Adding the foldchange & relative abudance 
Diff_genus_table_MD_T1D <- merge(Diff_genus_table_MD_T1D, DAA_genus, by = ".taxa")

#making the differential abundance tables for the article
genus_T1D_MD_article <- Diff_genus_table_MD_T1D[,-c(4,5,6,10, 11, 12,15)]
names(genus_T1D_MD_article)[7] <- "Log2FC"
names(genus_T1D_MD_article)[8] <- "Pvalue"
names(genus_T1D_MD_article)[9] <- "FDR"

#subsetting only the significant
genus_T1D_MD_article <- genus_T1D_MD_article %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)
View(genus_T1D_MD_article)

#Adding the foldchange & relative abundance
Diff_species_table_MD_T1D <- merge(Diff_species_table_MD_T1D, DAA_species, by = ".taxa")

#making the differential abundance tables for the article
Species_T1D_MD_article <- Diff_species_table_MD_T1D[,-c(4,5,6,10, 13, 14,15)]
names(Species_T1D_MD_article)[7] <- "Log2FC"
names(Species_T1D_MD_article)[8] <- "Pvalue"
names(Species_T1D_MD_article)[9] <- "FDR"

#subsetting only the significant
Species_T1D_MD_article <- Species_T1D_MD_article %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)
View(Species_T1D_MD_article)

#making boxplots of interesting taxa
taxa_table_genus <- taxa_table(MIDLOCMicrobiome_biomdata, rank = -2, unc = "drop", taxa = 0.000000001)
Sphingomonas_relabund <- subset(taxa_table_genus, .taxa %in% c("Bacteroides"))
Sphingomonas_relabund$percentage <- Sphingomonas_relabund$.abundance*100


Sphingomonas_plot <-ggplot(Sphingomonas_relabund, aes(x = group, y = percentage, color = group))+
  geom_violin(fill = NA, linewidth = 1.2)+
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.8)+
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.75)+
  theme_minimal()+
  labs(y = "Relative Abundance (%)")+
  ggtitle("Sphingomonas")
Sphingomonas_plot <- Sphingomonas_plot + theme(legend.position = "none")
Sphingomonas_plot

taxa_stats(MIDLOCMicrobiome_biomdata, taxa = "Sphingomonas", stat.by = "group", rank = -2, test = "wilcox", unc = "asis")

#Microbacterium
Microbacterium_relabund <- subset(taxa_table_genus, .taxa %in% c("Microbacterium"))
Microbacterium_relabund$percentage <- Microbacterium_relabund$.abundance*100


Microbacterium_plot <-ggplot(Microbacterium_relabund, aes(x = group, y = percentage, color = group))+
  geom_violin(fill = NA, linewidth = 1.2)+
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.8)+
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.75)+
  theme_minimal()+
  labs(y = "Relative Abundance (%)")+
  ggtitle("Microbacterium")
Microbacterium_plot <- Microbacterium_plot + theme(legend.position = "none")
Microbacterium_plot

#GGB9597
GGB9597_relabund <- subset(taxa_table_genus, .taxa %in% c("GGB9597"))
GGB9597_relabund$percentage <- GGB9597_relabund$.abundance*100


GGB9597_plot <-ggplot(GGB9597_relabund, aes(x = group, y = percentage, color = group))+
  geom_violin(fill = NA, linewidth = 1.2)+
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.8)+
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.75)+
  theme_minimal()+
  labs(y = "Relative Abundance (%)")+
  ggtitle("GGB9597")
GGB9597_plot <- GGB9597_plot + theme(legend.position = "none")
GGB9597_plot

#extra species plots
#c("Microbacterium_sp_T32", "Enterocloster_aldenensis", "GGB3677_SGB4990", "Clostridium_sp_Marseille_P3244", "GGB3109_SGB4121", "GGB9301_SGB14262")
species_relabund <- subset(species_table, .taxa %in% c("Bacteroides_uniformis"))
species_relabund$percentage <- species_relabund$.abundance*100


species_plot <-ggplot(species_relabund, aes(x = group, y = percentage, color = group))+
  geom_violin(fill = NA, linewidth = 1.2)+
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.8)+
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.75)+
  theme_minimal()+
  labs(y = "Relative Abundance (%)")+
  ggtitle("GGB9301_SGB14262")
species_plot <- species_plot + theme(legend.position = "none")
species_plot

taxa_stats(MIDLOCMicrobiome_biomdata, taxa = "Sphingomonas", stat.by = "group", rank = -2, test = "wilcox", unc = "asis")