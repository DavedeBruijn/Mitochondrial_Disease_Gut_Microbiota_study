#loading r packages
library(rbiom)
library(phyloseq)
library(vegan)
library(dplyr)

# setting library paths
MIDLOCmicrobioom_biompath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/04_cleaned/04_cleaned_MIDLOC.biom"

# import biom file into R
MIDLOCMicrobioom_biomdata <- as_rbiom(MIDLOCmicrobioom_biompath)
MIDLOCMicrobioom_biomdata

#create an otu_table 
otu_table <- as.matrix(MIDLOCMicrobioom_biomdata$counts)

#create an taxmap
taxmap <- taxa_map(MIDLOCMicrobioom_biomdata, unc = "asis")

#count the number of different taxa per rank
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", ".otu")

rank_counts_diff <- sapply(ranks, function(rank) {
  vals <- taxmap[[rank]]
  vals <- vals[vals != ""]
  length(unique(vals))
})

rank_counts_diff

#count the number of different taxa abundant in our data
nnz <- function (x) sum(x > 0)
percentage_rank_diff <- sapply(ranks, function(rank){
  (sample_apply(MIDLOCMicrobioom_biomdata, nnz, rank)/rank_counts_diff[rank])*100
})
mean <- colMeans(percentage_rank_diff)
percentage_rank_diff <- rbind(percentage_rank_diff, mean)

#count the number of taxa per rank with unassigned
rank_counts_all <- sapply(ranks, function(rank) {
  vals <- taxmap[[rank]]
  vals <- vals[vals != ""]
  length(vals)
})

rank_counts_all

#count the number of taxa abundant in our data
nnz <- function (x) sum(x > 0)
percentage_rank_all <- sapply(ranks, function(rank){
  (sample_apply(MIDLOCMicrobioom_biomdata, nnz, rank))})

percentage_rank_all

#total percentage per rank
rankstat <- sapply(taxmap[ranks], function(x) sum(x !=""))
percentage <- (rankstat/5879)*100
rankstat <- rbind(rankstat, percentage)

#taxa_stack
taxa_stacked(MIDLOCMicrobioom_biomdata, rank = "Phylum", order.by = ".sample", facet.by = "type")

#taxa_stats
taxa_table <- taxa_table(MIDLOCMicrobioom_biomdata, rank = "phylum", transform = "percent")
taxa_stats(MIDLOCMicrobioom_biomdata, rank = 2, stat.by = "type", test = "wilcox", transform = "percent")
taxa_boxplot(MIDLOCMicrobioom_biomdata, rank = "phylum", stat.by = "type", transform = "percent", y.transform = "none", taxa = 6)
taxa_boxplot(MIDLOCMicrobioom_biomdata, rank = "class", stat.by = "type", unc = "drop", transform = "percent")
taxa_boxplot(MIDLOCMicrobioom_biomdata, rank = "genus", stat.by = "type", unc = "drop", transform = "percent")

#Fimicutes and bacteroides table
Bact_relabund <- subset(taxa_table, .taxa %in% c("Bacteroidota"))
Firm_relabund <- subset(taxa_table, .taxa %in% c("Firmicutes"))
Bact_relabund$.abundance <- Bact_relabund$.abundance*100
Firm_relabund$.abundance <- Firm_relabund$.abundance*100

#Divided by timepoint
Bact_microbiom_MD <- subset(Bact_relabund, Bact_relabund$type %in% c("MD"))
Bact_microbiom_C <- subset(Bact_relabund, Bact_relabund$type %in% c("Control")) 
Bact_microbiom_T1D <- subset(Bact_relabund, Bact_relabund$type %in% c("T1D")) 
Firm_microbiom_MD <- subset(Firm_relabund, Firm_relabund$type %in% c("MD"))
Firm_microbiom_C <- subset(Firm_relabund, Firm_relabund$type %in% c("Control"))  
Firm_microbiom_T1D <- subset(Firm_relabund, Firm_relabund$type %in% c("T1D"))  

Firmicutes_diff_MD_C <- wilcox.test(Firm_microbiom_MD$.abundance, Firm_microbiom_C$.abundance, paired = FALSE)
Bacteroidota_diff_MD_C <- wilcox.test(Bact_microbiom_MD$.abundance , Bact_microbiom_C$.abundance, paired = FALSE)
Firmicutes_diff_MD_C
Bacteroidota_diff_MD_C

Firmicutes_diff_MD_T1D <- wilcox.test(Firm_microbiom_MD$.abundance, Firm_microbiom_T1D$.abundance, paired = FALSE)
Bacteroidota_diff_MD_T1D <- wilcox.test(Bact_microbiom_MD$.abundance , Bact_microbiom_T1D$.abundance, paired = FALSE)
Firmicutes_diff_MD_T1D
Bacteroidota_diff_MD_T1D

#plotting the firmicutes and bacteroides differece
library(ggplot2)
bacteroidota_plot <-ggplot(Bact_relabund, aes(x = type, y = .abundance, color = type))+
  geom_violin(fill = NA, linewidth = 1.2)+
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.8)+
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.75)+
  scale_color_manual(values = c("Control" = "#00BFC4", "MD" = "#F8766D", "T1D" = "#BA72FF"))+
  theme_minimal()+
  labs(y = "Relative Abundance (%)")+
  ggtitle("Bacteroidetes")+
  annotate("segment", x = 1, xend = 1.95, y = 83, yend = 83, linewidth = 0.5, color = "black") +
  annotate("segment", x = 1, xend = 1, y = 82, yend = 84, linewidth = 0.5, color = "black") +
  annotate("segment", x = 1.95, xend = 1.95, y = 82, yend = 84, linewidth = 0.5, color = "black") +
  annotate("text", x = 1.5, y = 84.75, label = "FDR = 0.0376", size = 2.5)+
  annotate("segment", x = 2.05, xend = 3, y = 83, yend = 83, linewidth = 0.5, color = "black") +
  annotate("segment", x = 2.05, xend = 2.05, y = 82, yend = 84, linewidth = 0.5, color = "black") +
  annotate("segment", x = 3, xend = 3, y = 82, yend = 84, linewidth = 0.5, color = "black") +
  annotate("text", x = 2.5, y = 84.75, label = "FDR = 0.0376", size = 2.5)
bacteroidota_plot <- bacteroidota_plot + theme(legend.position = "none")
bacteroidota_plot

Firmicutes_plot <-ggplot(Firm_relabund, aes(x = type, y = .abundance, color = type))+
  geom_violin(fill = NA, linewidth = 1.2)+
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.8)+
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.75)+
  scale_color_manual(values = c("Control" = "#00BFC4", "MD" = "#F8766D", "T1D" = "#BA72FF"))+
  theme_minimal()+
  labs(y = "Relative Abundance (%)")+
  ggtitle("Firmicutes")+
  annotate("segment", x = 2.05, xend = 3, y = 73, yend = 73, linewidth = 0.5, color = "black") +
  annotate("segment", x = 2.05, xend = 2.05, y = 72, yend = 74, linewidth = 0.5, color = "black") +
  annotate("segment", x = 3, xend = 3, y = 72, yend = 74, linewidth = 0.5, color = "black") +
  annotate("text", x = 2.5, y = 74.75, label = "FDR = 0.0320", size = 2.5)
Firmicutes_plot

library(patchwork)
Bact_Firm_plot <- bacteroidota_plot + Firmicutes_plot
Bact_Firm_plot + plot_annotation(
  tag_levels = "A",
  title = "Differences on Phylum rank between patients with MD and healthy controls",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

FDR_plot <- taxa_stats(MIDLOCMicrobioom_biomdata, rank = 2, stat.by = "type", test = "wilcox", taxa = c("Bacteroidota", "Firmicutes"), p.adj = "fdr")
FDR_plot
