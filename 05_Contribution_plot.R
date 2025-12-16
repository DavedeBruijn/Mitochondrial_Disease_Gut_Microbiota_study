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

#removing the extra rank
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

## species contribution towards the data.
## ASPASN-PWY: superpathway of L-aspartate and L-asparagine biosynthesis
## PWY-5103: L-isoleucine biosynthesis III
## BRANCHED-CHAIN-AA-SYN-PWY: superpathway of branched chain amino acid biosynthesis
## ILEUSYN-PWY: L-isoleucine biosynthesis I \\(from threonine\\)
## PWY-2942: L-lysine biosynthesis III
## PANTO-PWY: phosphopantothenate biosynthesis I
## PWY-7977: L-methionine biosynthesis IV
## SER-GLYSYN-PWY: superpathway of L-serine and glycine biosynthesis I
## THRESYN-PWY: superpathway of L-threonine biosynthesis
## PWY0-1586: peptidoglycan maturation \\(meso-diaminopimelate containing\\)
## PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I \\(bacteria\\)
## PYRIDNUCSYN-PWY: NAD de novo biosynthesis I \\(from aspartate\\)
## PWY-724: superpathway of L-lysine, L-threonine and L-methionine biosynthesis II
BCAA_otu_table <- otu_table[grep("^ASPASN-PWY: superpathway of L-aspartate and L-asparagine biosynthesis", rownames(otu_table)), ]
MIDLOC_pathway <- as_rbiom(BCAA_otu_table)
BCAA_otu_table <- as.data.frame(BCAA_otu_table)
view(BCAA_otu_table)

MIDLOC_pathway$metadata <- MIDLOCMicrobiome_biomdata$metadata
MIDLOC_pathway$taxonomy <- MIDLOCMicrobiome_biomdata$taxonomy
Diff_pathway_BCAA <- taxa_stats(MIDLOC_pathway, rank = 0, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop", transform = "none")
View(Diff_pathway_BCAA)

 # Reshape to long format for ggplot
BCAA_long <- BCAA_otu_table %>%
  tibble::rownames_to_column("species") %>%  
  pivot_longer(cols = -species,             
               names_to = "samples",  
               values_to = "value") 

#Add the groups to the long format
BCAA_long <- BCAA_long %>%
  mutate(
    # Extract sample number from the 'samples' column (assuming it has the format Sample1, Sample2, etc.)
    sample_num = as.numeric(gsub("MIDLOC", "", samples)),
    
    # Create the 'type' column based on the sample number
    type = case_when(
      sample_num >= 1 & sample_num <= 32 ~ "MD",       # Sample 1 to 32 are MD
      sample_num >= 33 & sample_num <= 92 ~ "Control",  # Sample 33 to 92 are Control
      sample_num >= 93 & sample_num <= 152 ~ "T1D",     # Sample 93 to 152 are T1D
    )
  )

BCAA_long$sample_num <- NULL

BCAA_long_MD_control <- subset(BCAA_long, BCAA_long$type == "MD" | BCAA_long$type == "Control")

results <- BCAA_long_MD_control %>%
  group_by(species) %>%
  summarise(
    wilcoxon_result = list(wilcox.test(value ~ type, data = cur_data())),
    mean_MD = mean(value[type == 'MD']),
    mean_Control = mean(value[type == 'Control']),
    Diff_C_MD = mean_Control - mean_MD,  # Difference between group means
    .groups = 'drop'
  )

results <- results %>%
  mutate(
    p_value = sapply(wilcoxon_result, function(x) x$p.value)
  )
results$wilcoxon_result <- NULL
results$FDR <- p.adjust(results$p_value , method = "fdr")
view(results)

#removing the pathway name so only the genus and species name remain
BCAA_long$species <- sub("^[^|]*\\|", "", BCAA_long$species)

# Aggregate data by group and species
BCAA_aggregated <- BCAA_long %>%
  group_by(type, species) %>%
  summarise(total_value = mean(value))

#plotting the contributions for the species
BCAA_aggregated <- BCAA_aggregated[-c(1,62,123),]

ggplot(BCAA_aggregated, aes(x = type, y = total_value, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Species Contribution to ASPASN-PWY: superpathway of L-aspartate and L-asparagine biosynthesis",
    x = "Groups",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "lines"),
        axis.title = element_text(size = 7),                  # Smaller axis title
        plot.title = element_text(size = 7, hjust = 0.0),     # Center title and adjust size
        legend.position = "right",                            # Position legend at the bottom
        legend.title = element_text(size = 5),                # Adjust legend title size
        legend.text = element_text(size = 4.8))

##
##species contribution towards the data. (PWY0-1586: peptidoglycan maturation \\(meso-diaminopimelate containing\\))
BCAA_otu_table <- otu_table[grep("^PWY0-1586: peptidoglycan maturation \\(meso-diaminopimelate containing\\)", rownames(otu_table)), ]
MIDLOC_pathway <- as_rbiom(BCAA_otu_table)
BCAA_otu_table <- as.data.frame(BCAA_otu_table)
view(BCAA_otu_table)

MIDLOC_pathway$metadata <- MIDLOCMicrobiome_biomdata$metadata
MIDLOC_pathway$taxonomy <- MIDLOCMicrobiome_biomdata$taxonomy
Diff_pathway_BCAA <- taxa_stats(MIDLOC_pathway, rank = 0, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop", transform = "none")
View(Diff_pathway_BCAA)

# Reshape to long format for ggplot
BCAA_long <- BCAA_otu_table %>%
  tibble::rownames_to_column("species") %>%  
  pivot_longer(cols = -species,             
               names_to = "samples",  
               values_to = "value") 

#Add the groups to the long format
BCAA_long <- BCAA_long %>%
  mutate(
    # Extract sample number from the 'samples' column (assuming it has the format Sample1, Sample2, etc.)
    sample_num = as.numeric(gsub("MIDLOC", "", samples)),
    
    # Create the 'type' column based on the sample number
    type = case_when(
      sample_num >= 1 & sample_num <= 32 ~ "MD",       # Sample 1 to 32 are MD
      sample_num >= 33 & sample_num <= 92 ~ "Control",  # Sample 33 to 92 are Control
      sample_num >= 93 & sample_num <= 152 ~ "T1D",     # Sample 93 to 152 are T1D
    )
  )

BCAA_long$sample_num <- NULL

BCAA_long_MD_control <- subset(BCAA_long, BCAA_long$type == "MD" | BCAA_long$type == "Control")

results <- BCAA_long_MD_control %>%
  group_by(species) %>%
  summarise(
    wilcoxon_result = list(wilcox.test(value ~ type, data = cur_data())),
    mean_MD = mean(value[type == 'MD']),
    mean_Control = mean(value[type == 'Control']),
    Diff_C_MD = mean_Control - mean_MD,  # Difference between group means
    .groups = 'drop'
  )

results <- results %>%
  mutate(
    p_value = sapply(wilcoxon_result, function(x) x$p.value)
  )
results$wilcoxon_result <- NULL
results$FDR <- p.adjust(results$p_value , method = "fdr")
view(results)

#removing the pathway name so only the genus and species name remain
BCAA_long$species <- sub("^[^|]*\\|", "", BCAA_long$species)

# Aggregate data by group and species
BCAA_aggregated <- BCAA_long %>%
  group_by(type, species) %>%
  summarise(total_value = mean(value))

#plotting the contributions for the species
BCAA_aggregated <- BCAA_aggregated[-c(1,86,171),]

ggplot(BCAA_aggregated, aes(x = type, y = total_value, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Species Contribution to PWY0-1586: peptidoglycan maturation (meso-diaminopimelate containing)",
    x = "Sample",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "lines"),
        axis.title = element_text(size = 7),                  # Smaller axis title
        plot.title = element_text(size = 7, hjust = 0.0),     # Center title and adjust size
        legend.position = "right",                            # Position legend at the bottom
        legend.title = element_text(size = 5),                # Adjust legend title size
        legend.text = element_text(size = 4.8))
##
##species contribution towards the data. (ASPASN-PWY: superpathway of L-aspartate and L-asparagine biosynthesis)
BCAA_otu_table <- otu_table[grep("^ASPASN-PWY: superpathway of L-aspartate and L-asparagine biosynthesis", rownames(otu_table)), ]
MIDLOC_pathway <- as_rbiom(BCAA_otu_table)
BCAA_otu_table <- as.data.frame(BCAA_otu_table)
view(BCAA_otu_table)

MIDLOC_pathway$metadata <- MIDLOCMicrobiome_biomdata$metadata
MIDLOC_pathway$taxonomy <- MIDLOCMicrobiome_biomdata$taxonomy
Diff_pathway_BCAA <- taxa_stats(MIDLOC_pathway, rank = 0, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop", transform = "none")
View(Diff_pathway_BCAA)

# Reshape to long format for ggplot
BCAA_long <- BCAA_otu_table %>%
  tibble::rownames_to_column("species") %>%  
  pivot_longer(cols = -species,             
               names_to = "samples",  
               values_to = "value") 

#Add the groups to the long format
BCAA_long <- BCAA_long %>%
  mutate(
    # Extract sample number from the 'samples' column (assuming it has the format Sample1, Sample2, etc.)
    sample_num = as.numeric(gsub("MIDLOC", "", samples)),
    
    # Create the 'type' column based on the sample number
    type = case_when(
      sample_num >= 1 & sample_num <= 32 ~ "MD",       # Sample 1 to 32 are MD
      sample_num >= 33 & sample_num <= 92 ~ "Control",  # Sample 33 to 92 are Control
      sample_num >= 93 & sample_num <= 152 ~ "T1D",     # Sample 93 to 152 are T1D
    )
  )

BCAA_long$sample_num <- NULL

BCAA_long_MD_control <- subset(BCAA_long, BCAA_long$type == "MD" | BCAA_long$type == "Control")

results <- BCAA_long_MD_control %>%
  group_by(species) %>%
  summarise(
    wilcoxon_result = list(wilcox.test(value ~ type, data = cur_data())),
    mean_MD = mean(value[type == 'MD']),
    mean_Control = mean(value[type == 'Control']),
    Diff_C_MD = mean_Control - mean_MD,  # Difference between group means
    .groups = 'drop'
  )

results <- results %>%
  mutate(
    p_value = sapply(wilcoxon_result, function(x) x$p.value)
  )
results$wilcoxon_result <- NULL
results$FDR <- p.adjust(results$p_value , method = "fdr")
view(results)

#removing the pathway name so only the genus and species name remain
BCAA_long$species <- sub("^[^|]*\\|", "", BCAA_long$species)

# Aggregate data by group and species
BCAA_aggregated <- BCAA_long %>%
  group_by(type, species) %>%
  summarise(total_value = mean(value))

#plotting the contributions for the species
BCAA_aggregated <- BCAA_aggregated[-c(1,62,123),]

ggplot(BCAA_aggregated, aes(x = type, y = total_value, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Species Contribution to ASPASN-PWY: superpathway of L-aspartate and L-asparagine biosynthesis",
    x = "Sample",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "lines"),
        axis.title = element_text(size = 7),                  # Smaller axis title
        plot.title = element_text(size = 7, hjust = 0.0),     # Center title and adjust size
        legend.position = "right",                            # Position legend at the bottom
        legend.title = element_text(size = 5),                # Adjust legend title size
        legend.text = element_text(size = 4.8))
