#Loading packages
library(rbiom)
library(phyloseq)
library(vegan)
library(dplyr)

# setting library paths
MIDLOCMicrobiome_biompath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/04_Cleaned/04_cleaned_MIDLOC.biom"

# import biom file into R
MIDLOC_microbiomedata <- as_rbiom(MIDLOCMicrobiome_biompath)
MIDLOC_microbiomedata

#removing the extra rank
colnames(MIDLOC_microbiomedata$taxonomy) <- c(".otu", "pathway", "genus", "species")
MIDLOC_microbiomedata

#normalize by the root since the colsum is not comparable
library(slam)
feature_table <- MIDLOC_microbiomedata$counts
col_sums <- col_sums(feature_table)
features_table_norm <- feature_table
features_table_norm$v <- features_table_norm$v/col_sums[feature_table$j]
col_sums(features_table_norm)
MIDLOC_microbiomedata$counts <- features_table_norm

#check the otu_table 
otu_table <- as.matrix(MIDLOC_microbiomedata$counts)
colSums(otu_table)
View(otu_table)


#Creating a PCoA plot with bray-Curtis
MIDLOCMicrobiome_PCoaplot<- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "peta", unc = "drop", rank = 0)
MIDLOCMicrobiome_PCoaplot

#Adding the centroids to the plot
ord_results <- MIDLOCMicrobiome_PCoaplot$data
centroids <- ord_results %>%
  group_by(MIDLOC_microbiomedata$metadata$group) %>%
  summarise(
    PC1 = mean(.x),
    PC2 = mean(.y))

#PCoA with vegan
dist_matrix <- vegdist(t(otu_table), method = "bray")
pcoa_results_vegdist <- cmdscale(dist_matrix, k=3)
pcoa_results_vegdist_eigenvalue <- cmdscale(dist_matrix, k=2, eig = TRUE)
pcoa_df_vegdist <- data.frame(PCoA1 = pcoa_results_vegdist[,1], PCoA2 = pcoa_results_vegdist[,2], PCoA3 = pcoa_results_vegdist[,3], eigenvalues = pcoa_results_vegdist_eigenvalue$eig, SampleID = rownames(pcoa_results_vegdist), SampleGroup = sample_data(MIDLOC_microbiomedata$metadata))

#eigenvalues and percentage of variance
eigenvalues <- pcoa_df_vegdist$eigenvalues
total_variance <- sum(eigenvalues)
pcoa_variance <- as.data.frame(eigenvalues, row.names = paste0("PC", 1:150))
pcoa_variance$variance <- pcoa_variance$eigenvalues/total_variance
pcoa_variance$percentage <- pcoa_variance$variance*100
pcoa_variance$percentage <- lapply(pcoa_variance$percentage, function(x) sprintf("%.1f", x))
View(pcoa_variance)

#making a nicer visual plot with ggplot
library(cowplot)
library(ggplot2)
colnames(centroids)[1] <- "group"
colnames(pcoa_df_vegdist)[7] <- "group"

MIDLOCMicrobiome_PCoaplot_show <-  bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "peta",unc = "drop", rank = 1) +
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = group),
             shape = 1, size = 1, stroke = 2) +
  xlab( paste("PC1 Explained dissimilarity", pcoa_variance$percentage[1] , "%")) + ylab( paste("PC2 Explained dissimilarity", pcoa_variance$percentage[2] , "%"))+
  theme(axis.ticks = element_line(color = "black"),axis.text = element_text(color = "black")) + 
  geom_hline(yintercept = 0, color = "grey") + geom_vline(xintercept = 0, color = "grey")+
  ggtitle("Gut Microbiome of patients with MD compared to healthy controls(Microbial)", subtitle = "weighted Bray-Curtis PCoA plot")
MIDLOCMicrobiome_PCoaplot_show

#adonis (permanova) 
PERMANOVA_results <- adonis2(dist_matrix~ pcoa_df_vegdist$group, by = 'margin')
PERMANOVA_results

#Pairwise Adonis2 analysis
library("OmicFlow")
pcoa_df_vegdist$group <- as.character(pcoa_df_vegdist$group)
pairwiseadonis_table <- pairwise_adonis(dist_matrix, pcoa_df_vegdist$group, perm = 9999)
pairwiseadonis_table

#Variance explained by each taxa
MIDLOCMicrobiome_PCoaplot_genus <-  bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "peta",unc = "drop", rank = 1, taxa = 6) +
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = group),
             shape = 1, size = 1, stroke = 2) +
  scale_color_manual(values = c("Control" = "#00BFC4", "MD" = "#F8766D", "T1D" = "#BA72FF"))+
  xlab( paste("PC1 Explained dissimilarity", pcoa_variance$percentage[1] , "%")) + ylab( paste("PC2 Explained dissimilarity", pcoa_variance$percentage[2] , "%"))+
  theme(axis.ticks = element_line(color = "black"),axis.text = element_text(color = "black")) + 
  geom_hline(yintercept = 0, color = "grey") + geom_vline(xintercept = 0, color = "grey")+
  ggtitle("Gut Microbiome of patients with MD compared to healthy controls(Microbial)", subtitle = "weighted Bray-Curtis PCoA plot")
MIDLOCMicrobiome_PCoaplot_genus

Taxavariancepca1 <- data.frame(MIDLOCMicrobiome_PCoaplot_genus$data$taxa_coords$.x, MIDLOCMicrobiome_PCoaplot_genus$data$taxa_coords$.label)
top_pc1 <- (Taxavariancepca1[,1]^2) / sum(Taxavariancepca1[,1]^2) * 100
Taxavariancepca1$top_pc1 <- top_pc1
Taxavariancepca1 <- Taxavariancepca1[order(-Taxavariancepca1$top_pc1),] 
print(Taxavariancepca1)

barplot( Taxavariancepca1$top_pc1, names.arg = Taxavariancepca1$ord.data.taxa_coords..label, main = "Genera dissimilarity of PC1", sub = "The genera who explain the most dissimilarity of PC1")

Taxavariancepca2 <- data.frame(MIDLOCMicrobiome_PCoaplot_genus$data$taxa_coords$.y, MIDLOCMicrobiome_PCoaplot_genus$data$taxa_coords$.label)
top_pc2 <- (Taxavariancepca2[,1]^2) / sum(Taxavariancepca2[,1]^2) * 100
Taxavariancepca2$top_pc2 <- top_pc2
Taxavariancepca2 <- Taxavariancepca2[order(-Taxavariancepca2$top_pc2),] 
print(Taxavariancepca2)

barplot( Taxavariancepca2$top_pc2, names.arg = Taxavariancepca2$ord.data.taxa_coords..label, main = "Genera dissimilarity of PC2", sub = "The genera who explain the most dissimilarity of PC2")

#Making a nicer plot with ggplot2
bdiv_microbial <- ggplot(pcoa_df_vegdist, aes(x = PCoA1, y = PCoA2, color = group))+
  geom_point(size = 1, alpha = 0.6)+
  theme_minimal()+
  geom_point(data = centroids, aes(x = PC1, y = PC2),
             shape = 18, size = 4)+
  scale_color_manual(values = c("Control" = "#00BFC4", "MD" = "#F8766D", "T1D" = "#BA72FF"))+
  stat_ellipse(type = "t", linetype = 1, alpha = 0.8)+
  xlab( paste("PC1,", pcoa_variance$percentage[1] , "% Dissimilarity")) + ylab( paste("PC2,", pcoa_variance$percentage[2] , "% Dissimilarity"))+
  theme(legend.direction = "horizontal")+
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16))+
  ggtitle("Weighted Bray-Curtis PCoA plot (Pathways)")
bdiv_microbial <- bdiv_microbial + theme(legend.position="bottom")
bdiv_microbial

#saving data to import figure data
library("readr")
centroids_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/05_Beta_diversity/centroids.tsv"
write_tsv(centroids, centroids_savepath)
pcoa_df_vegdist_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/05_Beta_diversity/pcoa_df_vegdist.tsv"
write_tsv(pcoa_df_vegdist, pcoa_df_vegdist_savepath)
pcoa_variance_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/05_Beta_diversity/pcoa_variance.tsv"
write_tsv(pcoa_variance ,pcoa_variance_savepath)
