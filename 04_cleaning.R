#Loading packages
library(rbiom)
library(phyloseq)
library(readr)
library(janitor)
library(tidyverse)
library(dplyr)

# setting library paths
MIDLOC_microbiompath <- 
MIDLOC_genefamiliespath <- 
MIDLOC_pathabundancepath <- 
MIDLOC_pathcoveragepath <- 
MIDLOC_metaphlan_tables <- 

# import biom file into R
MIDLOC_microbiomdata <- as_rbiom(MIDLOC_microbiompath)
MIDLOC_microbiomdata

#looking at the data types
glimpse(MIDLOC_microbiomdata$metadata)
glimpse(MIDLOC_microbiomdata$taxonomy)

#cleaning the taxonomy table by removing unnecessary characters (If you want to work with NA)
MIDLOC_microbiomdata$taxonomy <- MIDLOC_microbiomdata$taxonomy %>% 
  mutate(across(-1, ~str_remove(., "^.\\_\\_")))
View(MIDLOC_microbiomdata$taxonomy)

#filtering out the archea
unique(MIDLOC_microbiomdata$taxonomy$Kingdom)
features_table <- MIDLOC_microbiomdata$counts
class(features_table)

rows_to_keep <- MIDLOC_microbiomdata$taxonomy$Kingdom == "Bacteria"
features_table <- features_table[rows_to_keep,] 
MIDLOC_microbiomdata$counts <- features_table
unique(MIDLOC_microbiomdata$taxonomy$Kingdom)

#check the otu_table 
otu_table <- as.matrix(MIDLOC_microbiomdata$counts)
colSums(otu_table)
View(otu_table)

#normalize by the root since certain count are removed
library(slam)
col_sums <- col_sums(features_table)
features_table_norm <- features_table
features_table_norm$v <- features_table_norm$v*8/col_sums[features_table$j]
col_sums(features_table_norm)
MIDLOC_microbiomdata$counts <- features_table_norm

#renaming the column names
MIDLOC_microbiomdata$samples[31:90]<- paste0("MIDLOC", 33:92)
MIDLOC_microbiomdata$samples

withoutT0 <- sub("T0$", "", MIDLOC_microbiomdata$samples[1:30])
MIDLOC_microbiomdata$samples[1:30] <- withoutT0

#adding metadata 
MIDLOC_metadatapath <- 
#add metadata to R
metadata<-readr::read_tsv(MIDLOC_metadatapath)
print(metadata)

#could also rename if there are a lot of columns (sample is the name of the first column)
metadata <- rename(metadata, .sample = sample)

#add metadata to the biom file
MIDLOC_microbiomdata$metadata <- metadata
MIDLOC_microbiomdata

#cleaning with janitor
clean_metadata <- janitor::clean_names(MIDLOC_microbiomdata$metadata)
colnames(clean_metadata) <- gsub("^pat_", "", colnames(clean_metadata))
clean_metadata <- rename(clean_metadata, .sample = sample)
View(clean_metadata)

#adding the metadate back towards the biom file
MIDLOC_microbiomdata$metadata <- clean_metadata

#save rbiom object to a file
MIDLOC_biomdata_savepath <- 
write_biom(MIDLOC_microbiomdata, MIDLOC_biomdata_savepath)

#import tsv files
pathabundancetable <- data.table::fread(MIDLOC_pathabundancepath)
View(pathabundancetable)
pathcoveragetable <- data.table::fread(MIDLOC_pathcoveragepath)
View(pathcoveragetable)

#cleaning the first row of the TSV files
#Remove the number/letter code before the pathway name
cleaned_pathabundancetable <- pathabundancetable
cleaned_pathabundancetable$`# Pathway` <- sub("^.*?:\\s*", "",pathabundancetable$`# Pathway`)
cleaned_pathabundancetable$`# Pathway` <- sub("^\\s*", "",cleaned_pathabundancetable$`# Pathway`)
View(cleaned_pathabundancetable)

cleaned_pathcoveragetable <- pathcoveragetable
cleaned_pathcoveragetable$`# Pathway` <- sub("^.*?:\\s*", "",pathcoveragetable$`# Pathway` )
cleaned_pathcoveragetable$`# Pathway` <- sub("^\\s*", "",cleaned_pathcoveragetable$`# Pathway` )
View(cleaned_pathcoveragetable)

#clean the column names so there is only the sample_number
colnames(cleaned_genefamiliestable)[-1] <- sub("_.*", "",colnames(cleaned_genefamiliestable)[-1])
colnames(cleaned_pathabundancetable)[-1] <- sub("_.*", "",colnames(cleaned_pathabundancetable)[-1])
colnames(cleaned_pathcoveragetable)[-1] <- sub("_.*", "",colnames(cleaned_pathcoveragetable)[-1])

View(cleaned_genefamiliestable)
View(cleaned_pathabundancetable)
View(cleaned_pathcoveragetable)

#save tsv object to a file
pathabundancetable_savepath <- 
pathcoveragetable_savepath <- 

library(readr)
write_tsv(cleaned_pathabundancetable, pathabundancetable_savepath)
write_tsv(cleaned_pathcoveragetable, pathcoveragetable_savepath)


# check the distribution
hist(otu_table, ylim = range(0, 1000))
log_otu_table <- log10(otu_table + 1e-5)
hist(log_otu_table, ylim = range(0, 25000))