#Loading packages
library(rbiom)
library(phyloseq)
library(vegan)
library(dplyr)


# setting library paths
MIDLOCMicrobiome_biompath <- 

# import biom file into R
MIDLOCMicrobiome_biomdata <- as_rbiom(MIDLOCMicrobiome_biompath)
MIDLOCMicrobiome_biomdata

#create an otu_table 
otu_table <- as.matrix(MIDLOCMicrobiome_biomdata$counts)

#Differential abundance analysis
Diff_Family_table <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = -3, taxa = 0.0000001, stat.by = "type", test = "wilcox", unc = "drop")
Diff_genus_table <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = -2, taxa = 0.0000001, stat.by = "type", test = "wilcox", unc = "drop")
Diff_species_table <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = -1, taxa = 0.0000001, stat.by = "type", test = "wilcox", unc = "drop")
View(Diff_Family_table)
View(Diff_genus_table)
View(Diff_species_table)

#using OmicFlow to get the foldchange for the volcanoplot
library(OmicFlow)
library(data.table)

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

#using the foldchange
DAA_genus <- foldchange(data = genus_datatable, feature_rank = ".taxa", condition_A = c("_MD", "_C"), condition_B = c("_C", "_MD") ,paired = FALSE, condition_labels = names(genus_datatable)[-1])
View(DAA_genus)

#Adding the foldchange to the volcanoplot
Diff_genus_table <- merge(Diff_genus_table, DAA_genus, by = ".taxa")

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

#boxploting individual taxa
most_significantgenus <- top_hits$
rbiom::taxa_boxplot(MIDLOCMicrobiome_biomdata, rank = -2, taxa = most_significantgenus, stat.by = "type")
taxa_stats(MIDLOCMicrobiome_biomdata, rank = -2, taxa = most_significantgenus, stat.by = "type", test = "wilcox")
rbiom::taxa_boxplot(MIDLOCMicrobiome_biomdata, rank = -2, taxa = c("GGB9747", "Hydrogenoanaerobacterium"), stat.by = "type", transform = "percent")

Bifidobacterium <- c("Bifidobacterium")
taxa_boxplot(MIDLOCMicrobiome_biomdata, rank = -2, taxa = Bifidobacterium, stat.by = "type")
taxa_stats(MIDLOCMicrobiome_biomdata, rank = -2, taxa = Bifidobacterium, stat.by = "type" )

#oscillibacter
taxa_stats(MIDLOCMicrobiome_biomdata, rank = -2, taxa = "Oscillibacter", stat.by = "type")
taxa_stats(MIDLOCMicrobiome_biomdata, rank = -2, taxa = "Bacteroides", stat.by = "type")
taxa_stats(MIDLOCMicrobiome_biomdata, rank = -2, taxa = "Alistipes", stat.by = "type")

#making a species foldchange volcanoplot
#using OmicFlow to get the foldchange for the volcanoplot
#create a data.table 
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

#using the foldchange
DAA_species <- foldchange(data = species_datatable, feature_rank = ".taxa", condition_A = c("_MD", "_C"), condition_B = c("_C", "_MD") ,paired = FALSE, condition_labels = names(species_datatable)[-1])
View(DAA_species)

#Adding the foldchange to the volcanoplot
Diff_species_table <- merge(Diff_species_table, DAA_species, by = ".taxa")
View(Diff_species_table)

#calculating the relative abudance for the volcano plot
rel_abun_species <- species_table %>%
  group_by(.taxa) %>%
  summarise(rel_abun = mean(.abundance, na.rm = TRUE))
View(rel_abun_species)

#Adding the foldchange to the volcanoplot
Diff_species_table <- merge(Diff_species_table, rel_abun_species, by = ".taxa")

#making a volcano plot from the different abundance analysis Genus
library(ggplot2)
library(ggrepel)
volcanoplot_results_species <- Diff_species_table %>%
  mutate(log10_adj.pval = -log10(.adj.p),
         significant = .adj.p < 0.05)

top_hits_species <- volcanoplot_results_species %>%
  filter(.adj.p < 0.05) %>%
  arrange(.adj.p) 

hline_df <- data.frame(yintercept = -log10(0.05), Label = "p = 0.05")

DAA_species_plot <- ggplot(volcanoplot_results_species, aes(x = Log2FC_1, y = log10_adj.pval, colour = Log2FC_1))+
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
  geom_point(data = top_hits_species, aes(color = Log2FC_1, size = rel_abun), alpha = 1)+
  scale_size_continuous(range = c(0.1, 5), trans = "log10")+
  geom_text_repel(data = top_hits_species, aes(x = Log2FC_1, y = log10_adj.pval, label = .taxa), color = "black", size = 3, force = 2, max.overlaps = 20)
DAA_species_plot
library(patchwork)
DAA_species_plot + plot_annotation(
  title = "Difference in relative abundance of species in patients with MD and healthy controls",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)


