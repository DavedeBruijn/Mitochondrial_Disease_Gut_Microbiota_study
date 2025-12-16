#Loading packages
library(phyloseq)
library(readr)
library("dplyr")
library(rbiom)
library(tibble)

# setting library paths
MIDLOCMicrobiome_pathabundance_biompath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/01_loaded_data/01_loaded_MIDLOCpathabundancetable.tsv"

#import tsv files
pathabundancetable_Microbiome <- readr::read_tsv(MIDLOCMicrobiome_pathabundance_biompath)
View(pathabundancetable_Microbiome)

#taxonomy table
pathabundance_taxtable <- function(pathabundancetable_Microbiome) {
  pathway <- sub("\\|.*", "", pathabundancetable_Microbiome)
  
  genus <- sub(".*\\|g__([^\\.]+).*", "\\1", pathabundancetable_Microbiome)
  genus[!grepl("\\|", pathabundancetable_Microbiome)] <- ""
  
  species <- sub(".*\\.s__(.*)$", "\\1", pathabundancetable_Microbiome)
  species[!grepl("\\.s__", pathabundancetable_Microbiome)] <- ""
  
  tax_df <- data.frame(
    pathway = pathway,
    genus = genus,
    species = species,
    stringsAsFactors = FALSE
  )
  
  return(tax_df)
}
pathabundance_taxonomytable <- pathabundance_taxtable(pathabundancetable_Microbiome$`# Pathway`)
#convert to a matrix
pathabundancetable_Microbiome <- pathabundancetable_Microbiome %>% column_to_rownames(var = "# Pathway") %>% as.matrix()
View(pathabundancetable_Microbiome)

#convert to a biom file
MIDLOC_pathabundance_biom <- as_rbiom(pathabundancetable_Microbiome)

#check OTU table
otu_table <- as.matrix(MIDLOC_pathabundance_biom$counts)

#add taxonomy table
rownames(pathabundance_taxonomytable) <- rownames(pathabundancetable_Microbiome)
MIDLOC_pathabundance_biom$taxonomy <- pathabundance_taxonomytable
MIDLOC_pathabundance_biom

#Import metadata for subsetting
metadatasubsetting_path <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/02_subsetted_andBiom/metadata_subsetting.txt"
metadatasubsetting <- readr::read_tsv(metadatasubsetting_path)
metadatasubsetting <- rename(metadatasubsetting, .sample = sample)
metadatasubsetting$...2 <- NULL

#removing the _abundance so that the samples could be linked
without_Abundance <- sub("_Abundance$", "", MIDLOC_pathabundance_biom$samples)
MIDLOC_pathabundance_biom$samples <- without_Abundance

#adding the metadata
MIDLOC_pathabundance_biom$metadata <- metadatasubsetting

#subsetting to samples only used for MIDLOC diet analysis
MIDLOC_pathabundance_biom <- subset(MIDLOC_pathabundance_biom, MIDLOC_pathabundance_biom$metadata$Week %in% 0)
MIDLOC_pathabundance_biom
View(MIDLOC_pathabundance_biom$metadata)

#save rbiom object to a file
MIDLOCMicrobiome_biom_pathabundance_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/02_subsetted_andBiom/02_subsetted_andBiom_pathabundance.biom"
write_biom(MIDLOC_pathabundance_biom, MIDLOCMicrobiome_biom_pathabundance_savepath)

# setting library paths (Healthy controls)
MIDLOCMicrobiome_pathabundance_biompath_lifeline <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Lifelines data/01_loaded/01_loaded_MIDLOCpathabundancetable_lifeline.tsv"

#import tsv files
pathabundancetable_Microbiome_lifeline <- readr::read_tsv(MIDLOCMicrobiome_pathabundance_biompath_lifeline)
View(pathabundancetable_Microbiome_lifeline)

#taxonomy table
pathabundance_taxtable_lifeline <- function(pathabundancetable_Microbiome_lifeline) {
  pathway <- sub("\\|.*", "", pathabundancetable_Microbiome_lifeline)
  
  genus <- sub(".*\\|g__([^\\.]+).*", "\\1", pathabundancetable_Microbiome_lifeline)
  genus[!grepl("\\|", pathabundancetable_Microbiome_lifeline)] <- ""
  
  species <- sub(".*\\.s__(.*)$", "\\1", pathabundancetable_Microbiome_lifeline)
  species[!grepl("\\.s__", pathabundancetable_Microbiome_lifeline)] <- ""
  
  tax_df <- data.frame(
    pathway = pathway,
    genus = genus,
    species = species,
    stringsAsFactors = FALSE
  )
  
  return(tax_df)
}
pathabundance_taxonomytable_lifeline <- pathabundance_taxtable_lifeline(pathabundancetable_Microbiome_lifeline$`# Pathway`)
#convert to a matrix
pathabundancetable_Microbiome_lifeline <- pathabundancetable_Microbiome_lifeline %>% column_to_rownames(var = "# Pathway") %>% as.matrix()
View(pathabundancetable_Microbiome_lifeline)

#convert to a biom file
MIDLOC_pathabundance_biom_lifeline <- as_rbiom(pathabundancetable_Microbiome_lifeline)
MIDLOC_pathabundance_biom_lifeline
#check OTU table
otu_table <- as.matrix(MIDLOC_pathabundance_biom_lifeline$counts)

#add taxonomy table
rownames(pathabundance_taxonomytable_lifeline) <- rownames(otu_table)
MIDLOC_pathabundance_biom_lifeline$taxonomy <- pathabundance_taxonomytable_lifeline
MIDLOC_pathabundance_biom_lifeline

#save rbiom object to a file
MIDLOCMicrobiome_biom_pathabundance_lifeline_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/02_subsetted_andBiom/02_subsetted_andBiom_pathabundance_lifeline.biom"
write_biom(MIDLOC_pathabundance_biom_lifeline, MIDLOCMicrobiome_biom_pathabundance_lifeline_savepath)


# setting library paths (Type-1 diabetes)
MIDLOCMicrobiome_pathabundance_biompath_T1D <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/T1D data/01_loaded/01_loaded_MIDLOCpathabundancetable.tsv"

#import tsv files
pathabundancetable_Microbiome_T1D <- readr::read_tsv(MIDLOCMicrobiome_pathabundance_biompath_T1D)
View(pathabundancetable_Microbiome_T1D)

#only keeping the first 60 columns since the rest is double
pathabundancetable_Microbiome_T1D <- pathabundancetable_Microbiome_T1D[,1:61]

#taxonomy table
pathabundance_taxtable_T1D <- function(pathabundancetable_Microbiome_T1D) {
  pathway <- sub("\\|.*", "", pathabundancetable_Microbiome_T1D)
  
  genus <- sub(".*\\|g__([^\\.]+).*", "\\1", pathabundancetable_Microbiome_T1D)
  genus[!grepl("\\|", pathabundancetable_Microbiome_T1D)] <- ""
  
  species <- sub(".*\\.s__(.*)$", "\\1", pathabundancetable_Microbiome_T1D)
  species[!grepl("\\.s__", pathabundancetable_Microbiome_T1D)] <- ""
  
  tax_df <- data.frame(
    pathway = pathway,
    genus = genus,
    species = species,
    stringsAsFactors = FALSE
  )
  
  return(tax_df)
}
pathabundance_taxonomytable_T1D <- pathabundance_taxtable_T1D(pathabundancetable_Microbiome_T1D$`# Pathway`)
#convert to a matrix
pathabundancetable_Microbiome_T1D <- pathabundancetable_Microbiome_T1D %>% column_to_rownames(var = "# Pathway") %>% as.matrix()
View(pathabundancetable_Microbiome_T1D)

#convert to a biom file
pathabundancetable_Microbiome_T1D <- as_rbiom(pathabundancetable_Microbiome_T1D)
pathabundancetable_Microbiome_T1D
#check OTU table
otu_table <- as.matrix(pathabundancetable_Microbiome_T1D$counts)

#add taxonomy table
rownames(pathabundance_taxonomytable_T1D) <- rownames(otu_table)
pathabundancetable_Microbiome_T1D$taxonomy <- pathabundance_taxonomytable_T1D
pathabundancetable_Microbiome_T1D

#save rbiom object to a file
MIDLOCMicrobiome_biom_pathabundance_T1D_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/T1D data/02_metaphlan_subsetting/02_subsetted_andBiom_pathabundance_T1D.biom"
write_biom(pathabundancetable_Microbiome_T1D, MIDLOCMicrobiome_biom_pathabundance_T1D_savepath)

#extra checks for the pathcoverage and sanity check
filtered_pathabundance_MD <- pathabundancetable_Microbiome[!grepl("\\|", pathabundancetable_Microbiome$`# Pathway`),]
filtered_pathabundance_Control <- pathabundancetable_Microbiome_lifeline[!grepl("\\|", pathabundancetable_Microbiome_lifeline$`# Pathway`),]
filtered_pathabundance_T1D <- pathabundancetable_Microbiome_T1D[!grepl("\\|", pathabundancetable_Microbiome_T1D$`# Pathway`),]
