#Loading packages
library(rbiom)
library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)

# setting library paths
MIDLOCMicrobiome_biompath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/04_Cleaned/04_cleaned_MIDLOC.biom"

# import biom file into R
MIDLOCMicrobiome_biomdata <- as_rbiom(MIDLOCMicrobiome_biompath)
MIDLOCMicrobiome_biomdata

#renaming the rank
colnames(MIDLOCMicrobiome_biomdata$taxonomy) <- c(".otu", "pathway", "genus", "species")
MIDLOCMicrobiome_biomdata

#create an otu_table 
otu_table <- as.matrix(MIDLOCMicrobiome_biomdata$counts)
colSums(otu_table)

#normalize by the root since the colsum is not comparable
library(slam)
feature_table <- MIDLOCMicrobiome_biomdata$counts
col_sums <- col_sums(feature_table)
features_table_norm <- feature_table
features_table_norm$v <- features_table_norm$v/col_sums[feature_table$j]
col_sums(features_table_norm)
MIDLOCMicrobiome_biomdata$counts <- features_table_norm

#check the otu_table 
otu_table <- as.matrix(MIDLOCMicrobiome_biomdata$counts)
colSums(otu_table)
View(otu_table)

#Differential abundance analysis
Diff_pathway_table <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = 1, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_pathway_table)
Diff_pathway_table_genus <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = 2, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_pathway_table_genus)
Diff_pathway_table_species <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = 3, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_pathway_table_species)
Diff_pathway_table_otu <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = 0, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_pathway_table_otu)

#using OmicFlow to get the foldchange for the volcanoplot
library(OmicFlow)
library(data.table)

#create a data.table 
pathway_table <- taxa_table(MIDLOCMicrobiome_biomdata, rank = 1, taxa = 0.00000001, unc = "drop")
View(pathway_table)
pathway_list <- unique(pathway_table$.taxa)
View(pathway_list)

pathway_table <- as.data.table(pathway_table)
pathway_datatable <- dcast(pathway_table, .taxa ~ .sample, value.var = ".abundance")
View(pathway_datatable)

#Adding the label for the conditions
setnames(pathway_datatable, 
         old = names(pathway_datatable)[2:31], 
         new = paste0(names(pathway_datatable)[2:31], "_MD"))

# Rename columns 32 to 91 with _C suffix
setnames(pathway_datatable, 
         old = names(pathway_datatable)[32:91], 
         new = paste0(names(pathway_datatable)[32:91], "_C"))
# Rename columns 92 to 151 with _C suffix
setnames(pathway_datatable, 
         old = names(pathway_datatable)[92:151], 
         new = paste0(names(pathway_datatable)[92:151], "_T1D"))

#using the foldchange
DAA_pathway <- foldchange(data = pathway_datatable, feature_rank = ".taxa", condition_A = c("_MD", "_MD"), condition_B = c("_C", "_T1D"), paired = FALSE, condition_labels = names(pathway_datatable)[-1])
View(DAA_pathway)

#Adding the foldchange to the volcanoplot
Diff_pathway_table <- merge(Diff_pathway_table, DAA_pathway, by = ".taxa")

#calculating the relative abudance for the volcano plot
rel_abun_genus <- genus_table %>%
  group_by(.taxa) %>%
  summarise(rel_abun = mean(.abundance, na.rm = TRUE))
View(rel_abun_genus)

#Adding the foldchange to the volcanoplot
Diff_genus_table <- merge(Diff_genus_table, rel_abun_genus, by = ".taxa")

#making a volcano plot from the different abundance analysis Genus
library(ggplot2)
library(ggrepel)
volcanoplot_results <- Diff_genus_table %>%
  mutate(log10_adj.pval = -log10(.adj.p),
         significant = .adj.p < 0.05)

top_hits <- volcanoplot_results %>%
  filter(.adj.p < 0.05) %>%
  arrange(.adj.p) 

hline_df <- data.frame(yintercept = -log10(0.05), Label = "p = 0.05")

DAA_genus_plot <- ggplot(volcanoplot_results, aes(x = Log2FC_1, y = log10_adj.pval, colour = Log2FC_1))+
  geom_point(alpha = 1, size = 1.5) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
  theme_minimal()+
  labs(x = "Fold Change log2(prior/after)",
       y = "-log10(p-value)",
       size = "Relative Abundance",
       colour = "Fold Change")+
  geom_hline(data = hline_df, aes(yintercept = yintercept), color = "blue", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = -23, y = -log10(0.030), label = "p < 0.05", color = "black", size = 3)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = top_hits, aes(color = Log2FC_1, size = rel_abun), alpha = 1)+
  scale_size_continuous(range = c(0.1, 5), trans = "log10")+
  geom_text_repel(data = top_hits, aes(x = Log2FC_1, y = log10_adj.pval, label = .taxa), color = "black", size = 3, force = 2, max.overlaps = 20)
DAA_genus_plot
library(patchwork)
DAA_genus_plot + plot_annotation(
  title = "Difference in relative abundance of genus in patients with MD and healthy controls",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)
DAA_genus_plot