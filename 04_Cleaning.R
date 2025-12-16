#Loading packages
library(rbiom)
library(phyloseq)
library(readr)
library(janitor)
library(tidyverse)
library(dplyr)

# setting library paths
MIDLOC_pathabundancepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/03_Merged/03_Merged.biom"
MIDLOC_pathcoveragepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Lifelines data/read_annotation/merged_pathcoverage.tsv"

# import biom file into R
MIDLOC_microbiomdata_pathabundance <- as_rbiom(MIDLOC_pathabundancepath)
MIDLOC_microbiomdata_pathabundance

#renaming the removed name
colnames(MIDLOC_microbiomdata_pathabundance$taxonomy) <- c(".otu", "pathway", "genus", "species")
MIDLOC_microbiomdata_pathabundance

#looking at the data types
view(MIDLOC_microbiomdata_pathabundance$metadata)
View(MIDLOC_microbiomdata_pathabundance$taxonomy)

#check the otu_table 
otu_table <- as.matrix(MIDLOC_microbiomdata_pathabundance$counts)
colSums(otu_table)
View(otu_table)

#renaming the column names
MIDLOC_microbiomdata_pathabundance$samples[1:9]<- paste0("MIDLOC0", 10:18)
MIDLOC_microbiomdata_pathabundance$samples[10:16]<- paste0("MIDLOC00", 1:7)
MIDLOC_microbiomdata_pathabundance$samples[17]<- paste0("MIDLOC00", 9)
MIDLOC_microbiomdata_pathabundance$samples[18:28]<- paste0("MIDLOC0", 19:29)
MIDLOC_microbiomdata_pathabundance$samples[29:30]<- paste0("MIDLOC0", 31:32)
MIDLOC_microbiomdata_pathabundance$samples[31:97]<- paste0("MIDLOC0", 33:99)
MIDLOC_microbiomdata_pathabundance$samples[98:150]<- paste0("MIDLOC", 100:152)
MIDLOC_microbiomdata_pathabundance$samples

#adding metadata 
MIDLOC_metadatapath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/metdata.tsv"
#add metadata to R
metadata<-readr::read_tsv(MIDLOC_metadatapath)
print(metadata)

#could also rename if there are a lot of columns (sample is the name of the first column)
metadata <- rename(metadata, .sample = sample)

#add metadata to the biom file
MIDLOC_microbiomdata_pathabundance$metadata <- metadata
MIDLOC_microbiomdata_pathabundance

#cleaning with janitor
clean_metadata <- janitor::clean_names(MIDLOC_microbiomdata_pathabundance$metadata)
colnames(clean_metadata) <- gsub("^pat_", "", colnames(clean_metadata))
clean_metadata <- rename(clean_metadata, .sample = sample)
clean_metadata$laxantia_21 <- NULL
colnames(clean_metadata)[c(2,10,11,14)] <- c("group", "diabetes", "diabetes_code", "laxantia") 
View(clean_metadata)

#adding the metadate back towards the biom file
MIDLOC_microbiomdata_pathabundance$metadata <- clean_metadata

#save rbiom object to a file 
MIDLOC_biomdata_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/04_Cleaned/04_cleaned_MIDLOC.biom"
write_biom(MIDLOC_microbiomdata_pathabundance, MIDLOC_biomdata_savepath)

# check the distribution
hist(otu_table, ylim = range(0, 1000))
log_otu_table <- log10(otu_table + 1e-5)
hist(log_otu_table, ylim = range(0, 25000))