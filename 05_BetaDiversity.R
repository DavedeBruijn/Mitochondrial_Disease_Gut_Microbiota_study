#Loading packages
library(rbiom)
library(phyloseq)
library(vegan)
library(dplyr)

# setting library paths
MIDLOCMicrobiome_biompath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/04_cleaned/04_cleaned_MIDLOC.biom"

# import biom file into R
MIDLOC_microbiomedata <- as_rbiom(MIDLOCMicrobiome_biompath)
MIDLOC_microbiomedata

#create an otu_table 
otu_table <- as.matrix(MIDLOC_microbiomedata$counts)
View(otu_table)

names(MIDLOC_microbiomedata$metadata)[2] <- "group"

#normalize by the root since certain count were removed
library(slam)
features_table <- MIDLOC_microbiomedata$counts
col_sums <- col_sums(features_table)
features_table_norm <- features_table
features_table_norm$v <- features_table_norm$v/col_sums[features_table$j]
col_sums(features_table_norm)
MIDLOC_microbiomedata$counts <- features_table_norm

#Creating a PCoA plot with bray-Curtis
MIDLOCMicrobiome_PCoaplot<- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "peta", unc = "drop", rank = -2)
MIDLOCMicrobiome_PCoaplot

#Adding the centroids to the plot
ord_results <- MIDLOCMicrobiome_PCoaplot$data
centroids <- ord_results %>%
  group_by(MIDLOC_microbiomedata$metadata$group) %>%
  summarise(
    PC1 = mean(.x),
    PC2 = mean(.y))

#adding the centroids of all pcoa
centroids <- pcoa_df_vegdist %>%
  group_by(pcoa_df_vegdist$group) %>%
  summarise(
    PC1 = mean(PCoA1),
    PC2 = mean(PCoA2),
    PC3 = mean(PCoA3))


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

MIDLOCMicrobiome_PCoaplot_show <-  bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "peta",unc = "drop", rank = -2) +
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
MIDLOCMicrobiome_PCoaplot_genus <-  bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "peta",unc = "drop", rank = -2, taxa = 0.01) +
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = group),
             shape = 1, size = 1, stroke = 2) +
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
  geom_point(size = 1, alpha = 0.8)+
  theme_minimal()+
  geom_point(data = centroids, aes(x = PC1, y = PC2),
             shape = 18, size = 5)+
  scale_color_manual(values = c("Control" = "#00BFC4", "MD" = "#F8766D", "T1D" = "#BA72FF"))+
  stat_ellipse(type = "t", linetype = 1, alpha = 1)+
  xlab( paste("PC1,", pcoa_variance$percentage[1] , "% Dissimilarity")) + ylab( paste("PC2,", pcoa_variance$percentage[2] , "% Dissimilarity"))+
  theme(legend.direction = "horizontal")+
  theme(
    axis.title.x = element_text(size = 10, face = "bold"),  
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    title = element_text(size = 10))+
  ggtitle("Weighted Bray-Curtis PCoA plot (Microbial)")
bdiv_microbial <- bdiv_microbial + theme(legend.position="bottom")
bdiv_microbial

#using pcoa3 for the plot
#Making a nicer plot with ggplot2
bdiv_microbial_pcoa3 <- ggplot(pcoa_df_vegdist, aes(x = PCoA1, y = PCoA3, color = group))+
  geom_point(size = 1, alpha = 0.8)+
  theme_minimal()+
  geom_point(data = centroids, aes(x = PC1, y = PC3),
             shape = 18, size = 3)+
  scale_color_manual(values = c("Control" = "#00BFC4", "MD" = "#F8766D", "T1D" = "#BA72FF"))+
  stat_ellipse(type = "t", linetype = 1, alpha = 1)+
  xlab( paste("PC1,", pcoa_variance$percentage[1] , "% Dissimilarity")) + ylab( paste("PC3,", pcoa_variance$percentage[3] , "% Dissimilarity"))+
  theme(legend.direction = "horizontal")+
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16))+
  ggtitle("Weighted Bray-Curtis PCoA plot (Microbial)")
bdiv_microbial_pcoa3 <- bdiv_microbial_pcoa3 + theme(legend.position="bottom")
bdiv_microbial_pcoa3

#looking at the effect of fibre intake or other variables
#Creating a PCoA plot with bray-Curtis
MIDLOCMicrobiome_PCoaplot_fibre <- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "pem", unc = "drop", rank = -2)+ geom_text(aes(label = MIDLOC_microbiomedata$metadata$fibre))
MIDLOCMicrobiome_PCoaplot_fibre

MIDLOCMicrobiome_PCoaplot_NMDAS <- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "pem", unc = "drop", rank = -2)+ geom_text(aes(label = MIDLOC_microbiomedata$metadata$nmdas))
MIDLOCMicrobiome_PCoaplot_NMDAS

MIDLOCMicrobiome_PCoaplot_diabetes <- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "pem", unc = "drop", rank = -2)+ geom_text(aes(label = MIDLOC_microbiomedata$metadata$diabetes))
MIDLOCMicrobiome_PCoaplot_diabetes

MIDLOCMicrobiome_PCoaplot_laxantia <- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "pem", unc = "drop", rank = -2)+ geom_text(aes(label = MIDLOC_microbiomedata$metadata$laxantia))
MIDLOCMicrobiome_PCoaplot_laxantia

MIDLOCMicrobiome_PCoaplot_constipation <- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "pem", unc = "drop", rank = -2)+ geom_text(aes(label = MIDLOC_microbiomedata$metadata$constipation))
MIDLOCMicrobiome_PCoaplot_constipation

MIDLOCMicrobiome_PCoaplot_diarhee <- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "pem", unc = "drop", rank = -2)+ geom_text(aes(label = MIDLOC_microbiomedata$metadata$diarrhea))
MIDLOCMicrobiome_PCoaplot_diarhee

MIDLOCMicrobiome_PCoaplot_sex <- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "pem", unc = "drop", rank = -2)+ geom_text(aes(label = MIDLOC_microbiomedata$metadata$sex))
MIDLOCMicrobiome_PCoaplot_sex

MIDLOCMicrobiome_PCoaplot_BMI <- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "pem", unc = "drop", rank = -2)+ geom_text(aes(label = MIDLOC_microbiomedata$metadata$bmi))
MIDLOCMicrobiome_PCoaplot_BMI

MIDLOCMicrobiome_PCoaplot_probiotica <- bdiv_ord_plot(MIDLOC_microbiomedata, bdiv = "Bray-Curtis", ord = "PCoA", stat.by = "group", layers = "pem", unc = "drop", rank = -2)+ geom_text(aes(label = MIDLOC_microbiomedata$metadata$probiotica))
MIDLOCMicrobiome_PCoaplot_probiotica


#betadisper to look if variance could be a possible cofounder
bd <- betadisper(dist_matrix, pcoa_df_vegdist$group)
anova(bd)
permutest(bd)

#Combining the Beta Diversity plots for the article
#adding the functional plot for patchwork to combine them
centroids_functional <- data.table::fread("C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/05_Beta_diversity/centroids.tsv")
pcoa_df_vegdist_functional <- data.table::fread("C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/05_Beta_diversity/pcoa_df_vegdist.tsv")
pcoa_variance_functional <- data.table::fread("C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/05_Beta_diversity/pcoa_variance.tsv")
pcoa_variance_functional$percentage <- pcoa_variance_functional$variance*100
pcoa_variance_functional$percentage <- lapply(pcoa_variance_functional$percentage, function(x) sprintf("%.1f", x))

#Making the PCoA plot for the functional
bdiv_functional <- ggplot(pcoa_df_vegdist_functional, aes(x = PCoA1, y = PCoA2, color = group))+
  geom_point(size = 1, alpha = 0.8)+
  theme_minimal()+
  geom_point(data = centroids_functional, aes(x = PC1, y = PC2),
             shape = 18, size = 5)+
  scale_color_manual(values = c("Control" = "#00BFC4", "MD" = "#F8766D", "T1D" = "#BA72FF"))+
  stat_ellipse(type = "t", linetype = 1, alpha = 1)+
  xlab( paste("PC1,", pcoa_variance_functional$percentage[1] , "% Dissimilarity")) + ylab( paste("PC2,", pcoa_variance_functional$percentage[2] , "% Dissimilarity"))+
  theme(legend.direction = "horizontal")+
  theme(
    axis.title.x = element_text(size = 10, face = "bold"),  
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    title = element_text(size = 10))+
  ggtitle("Weighted Bray-Curtis PCoA plot (Metacyc Pathways)")
bdiv_functional <- bdiv_functional + theme(legend.position="bottom")
bdiv_functional

#adding the redundancy analysis plot
RDAresults <- data.table::fread("C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/05_RDA/RDAresults.tsv")

RDAresults <- RDAresults %>%
  mutate(Variable = factor(Variable, levels = Variable[order(r.squared, decreasing = FALSE)]))

RDA_plot <- ggplot(RDAresults, aes(x = ExplainedPercentage, y = Variable, fill = color))+
  geom_col(alpha =  1)+
  scale_fill_identity(name = "Legend",              
                      labels = c("Significant", "Not Significant"),   
                      breaks = c("brown1", "cornflowerblue"), 
                      guide = "legend")+
  labs(x = "Explained dissimilarity (%)", y = "Study variables")+
  theme_minimal()+
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"))+
  ggtitle("(Microbial)")
RDA_plot <- RDA_plot + theme(legend.position= "bottom")

#functionality
RDAresults_func <- data.table::fread("C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/05_RDA/RDAresults.tsv")

RDAresults_func <- RDAresults_func %>%
  mutate(Variable = factor(Variable, levels = Variable[order(r.squared, decreasing = FALSE)]))

RDA_plot_func <- ggplot(RDAresults_func, aes(x = ExplainedPercentage, y = Variable, fill = color))+
  geom_col(alpha =  1)+
  scale_fill_identity(name = "Legend",              
                      labels = c("Significant", "Not Significant"),   
                      breaks = c("brown1", "cornflowerblue"), 
                      guide = "legend")+
  labs(x = "Explained dissimilarity (%)", y = "Study variables")+
  theme_minimal()+
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"))+
  ggtitle("(Metacyc Pathways)")
RDA_plot_func <- RDA_plot_func + theme(legend.position= "bottom")


#combinding the plots with patchwork
library(patchwork)
bdiv_article <- bdiv_microbial + bdiv_functional+ RDA_plot + RDA_plot_func + plot_annotation(
  tag_levels = "A"
)
bdiv_article
