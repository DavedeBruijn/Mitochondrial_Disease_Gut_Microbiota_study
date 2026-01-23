#Loading packages
library(rbiom)
library(phyloseq)
library(vegan)
library(dplyr)
library(slam)
library(data.table)

# setting library paths
MIDLOCMicrobiome_biompath <- 

# import biom file into R
MIDLOCMicrobiome_biomdata <- as_rbiom(MIDLOCMicrobiome_biompath)
MIDLOCMicrobiome_biomdata

#renaming the ranks
colnames(MIDLOCMicrobiome_biomdata$taxonomy) <- c(".otu", "pathway", "genus", "species")
MIDLOCMicrobiome_biomdata

#create an otu_table 
otu_table <- as.matrix(MIDLOCMicrobiome_biomdata$counts)
colSums(otu_table)

#changing type to group
names(MIDLOCMicrobiome_biomdata$metadata)[2] <- "group"

#normalize by the root since certain count were removed
features_table <- MIDLOCMicrobiome_biomdata$counts
col_sums <- col_sums(features_table)
features_table_norm <- features_table
features_table_norm$v <- features_table_norm$v/col_sums[features_table$j]
col_sums(features_table_norm)
MIDLOCMicrobiome_biomdata$counts <- features_table_norm

#check new otu_table 
otu_table <- as.matrix(MIDLOCMicrobiome_biomdata$counts)
colSums(otu_table)

#adding the pathcoverage 30 table
pathcoverage_30_savepath <- 
pathcoverage_30 <- readr::read_tsv(pathcoverage_30_savepath)

otu_table_30 <- otu_table[row.names(otu_table) %in% pathcoverage_30$X..Pathway,]

#adjusting the otu_table to the biom file
MIDLOCMicrobiome_biomdata$counts <- otu_table_30

#Differential abundance analysis
Diff_pathway_table_genus <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = 2, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_pathway_table_genus)
Diff_pathway_table_species <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = 3, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_pathway_table_species)
Diff_pathway_table_otu <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = 0, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_pathway_table_otu)

#adding the filtered pathcoverage 50 table
pathcoverage_50_savepath <- 
pathcoverage_50 <- readr::read_tsv(pathcoverage_50_savepath)

otu_table_50 <- otu_table[row.names(otu_table) %in% pathcoverage_50$X..Pathway,]

#adjusting the otu_table to the biom file
MIDLOCMicrobiome_biomdata$counts <- otu_table_50

#Differential abundance analysi
Diff_pathway_table <- taxa_stats(MIDLOCMicrobiome_biomdata, rank = 1, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
colnames(Diff_pathway_table)[1] <- "pathways"
View(Diff_pathway_table)

#create a data.table 
pathway_table <- taxa_table(MIDLOCMicrobiome_biomdata, rank = 1, taxa = 0.00000001, unc = "drop")
View(pathway_table)
pathway_list <- unique(pathway_table$.taxa)
View(pathway_list)

pathway_table <- data.table::as.data.table(pathway_table)
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
library("OmicFlow")
DAA_pathway <- foldchange(data = pathway_datatable, feature_rank = ".taxa", condition_A = c("_MD", "_MD"), condition_B = c("_C", "_T1D") ,paired = FALSE, condition_labels = names(pathway_datatable)[-1])
View(DAA_pathway)

DAA_pathway$FDR_1 <-p.adjust(DAA_pathway$pvalue_1 , method = "fdr")
DAA_pathway$FDR_2 <-p.adjust(DAA_pathway$pvalue_2 , method = "fdr")

#checking with OmicFlow
pathcoverage_50_savepath <- 
pathcoverage_50 <- readr::read_tsv(pathcoverage_50_savepath)

otu_table_50 <- otu_table[row.names(otu_table) %in% pathcoverage_50$X..Pathway,]

#adjusting the otu_table to the biom file
MIDLOCMicrobiome_biomdata$counts <- otu_table_50
#changing type to group
metadata <- as.data.table(MIDLOCMicrobiome_biomdata$metadata)
countData <- as.matrix(MIDLOCMicrobiome_biomdata$counts)
featureData <- as.data.frame(MIDLOCMicrobiome_biomdata$taxonomy)

colnames(metadata)[1]<- "SAMPLE_ID"
colnames(featureData) <- c("FEATURE_ID","Pathway","Genus", "Species")
featureData$Kingdom <- NA
featureData$Phylum  <- NA
featureData$Class   <- NA
featureData$Order   <- NA
featureData$Family  <- NA


names(metadata)[2] <- "CONTRAST_Group"
metadata$CONTRAST_Group <- as.character(metadata$CONTRAST_Group)

taxa <- metagenomics$new(
  metaData = metadata,
  countData = countData,
  featureData = featureData
)
res <- taxa$DFE(
  feature_rank = "Pathway",
  paired = FALSE,
  normalize = FALSE,
  condition.group = "CONTRAST_Group",
  condition_A = "MD",
  condition_B = "Control"
)

View(res$data)

#calculating the relative abudance for the volcano plot
rel_abun <- as.data.frame(res$data)
colnames(DAA_pathway)[1] <- "Pathway"

DAA_pathway <- merge(DAA_pathway, rel_abun[,c("Pathway", "rel_abun")], by = "Pathway", all.x = TRUE)
DAA_pathway$rel_abun <- (DAA_pathway$rel_abun/sum(DAA_pathway$rel_abun))*100

#making a volcano plot from the different abundance analysis Genus
volcanoplot_results <- DAA_pathway %>%
  mutate(log10_adj.pval = -log10(FDR_1),
         significant = FDR_1 < 0.05)

top_hits <- volcanoplot_results %>%
  filter(FDR_1 < 0.05) %>%
  arrange(FDR_1) 

#removing the code for readability
top_hits$Pathway <- sub("^.*?:\\s*", "", top_hits$Pathway)

hline_df <- data.frame(yintercept = -log10(0.05), Label = "FDR = 0.05")

library(ggplot2)
library(ggrepel)
DAA_pathway_plot <- ggplot(volcanoplot_results, aes(x = Log2FC_1, y = log10_adj.pval, colour = Log2FC_1))+
  geom_point(alpha = 1, size = 2) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
  theme_minimal()+
  labs(x = "Fold Change log2(MD/Control)",
       y = "-log10(FDR)",
       size = "Relative Abundance (%)",
       colour = "Fold Change",
       title = "Pathway")+
  geom_hline(data = hline_df, aes(yintercept = yintercept), color = "blue", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = -0.5, y = -log10(0.035), label = "FDR < 0.05", color = "black", size = 3)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = top_hits, aes(color = Log2FC_1, size = rel_abun), alpha = 1)+
  scale_size_continuous(range = c(0.1, 4), trans = "log10",breaks = c(10, 1, 0.1, 0.01), labels = c("10", "1", "0.1", "0.01"))+
  geom_text_repel(data = top_hits, aes(x = Log2FC_1, y = log10_adj.pval, label = Pathway), color = "black", size = 2.5, force = 2.5, max.overlaps = 20)
DAA_pathway_plot

# Diagram image
library(cowplot)
diagram <- ggdraw() + draw_image()

p2 <- ggplot() + geom_blank() + theme_minimal()

# Combine
plot_grid(DAA_pathway_plot, p2, ncol = 1, labels = c("A", "B"))

#making the tables for the article
#making the Control - MD comparison
MIDLOC_C_MD <- subset(MIDLOCMicrobiome_biomdata, MIDLOCMicrobiome_biomdata$metadata$group != "T1D")
MIDLOC_C_MD

Diff_pathway_table_C_MD <- taxa_stats(MIDLOC_C_MD, rank = 1, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_pathway_table_C_MD)

#Adding the foldchange to tables
colnames(Diff_pathway_table_C_MD)[1] <- "Pathway"
Diff_pathway_table_C_MD <- merge(Diff_pathway_table_C_MD, DAA_pathway, by = "Pathway")

#making the differential abundance tables for the article
Pathway_C_MD_article <- Diff_pathway_table_C_MD[,-c(4,5,6,10, 13, 14,16)]
names(Pathway_C_MD_article)[7] <- "Log2FC"
names(Pathway_C_MD_article)[8] <- "Pvalue"
names(Pathway_C_MD_article)[9] <- "FDR"

#subsetting only the significant
Pathway_C_MD_article <- Pathway_C_MD_article %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)
View(Pathway_C_MD_article)

#
#making the Type 1 Diabetes - MD comparison
MIDLOC_MD_T1D <- subset(MIDLOCMicrobiome_biomdata, MIDLOCMicrobiome_biomdata$metadata$group != "Control")
MIDLOC_MD_T1D

#swap the levels of the groups so that it becomes: MD control
MIDLOC_MD_T1D$metadata$group <- factor(MIDLOC_MD_T1D$metadata$group, levels = rev(levels(MIDLOC_MD_T1D$metadata$group)))


#Differential abundance analysis table MD_T1D
Diff_pathway_table_MD_T1D <- taxa_stats(MIDLOC_MD_T1D, rank = 1, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_pathway_table_MD_T1D)


#Adding the foldchange & relative abudance 
colnames(Diff_pathway_table_MD_T1D)[1] <- "Pathway"
Diff_pathway_table_MD_T1D <- merge(Diff_pathway_table_MD_T1D, DAA_pathway, by = "Pathway")

#making the differential abundance tables for the article
pathway_T1D_MD_article <- Diff_pathway_table_MD_T1D[,-c(4,5,6,10, 11, 12,15)]
names(pathway_T1D_MD_article)[7] <- "Log2FC"
names(pathway_T1D_MD_article)[8] <- "Pvalue"
names(pathway_T1D_MD_article)[9] <- "FDR"

#subsetting only the significant
pathway_T1D_MD_article <- pathway_T1D_MD_article %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)
View(pathway_T1D_MD_article)

