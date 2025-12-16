#Loading packages
library("biomformat")
library(rbiom)

#adding metaphlan to the lifelines data
#setting library paths Lifeline data
MIDLOC_biompath_Lifeline <- 
MIDLOC_metaphlan_tables_Lifeline <- 

# import biom files into R
MIDLOC_biomdata_lifeline <- as_rbiom(MIDLOC_biompath_Lifeline)
MIDLOC_biomdata_lifeline

#adding otu names towards the otu's
#import tsv files
metaphlan_tables_lifeline <- data.table::fread(MIDLOC_metaphlan_tables_Lifeline)
View(metaphlan_tables_lifeline)

#change the otu_1 numbering towards the otu names
MIDLOC_biomdata_lifeline$otus <- metaphlan_tables_lifeline$V1
otu_table_lifeline <- as.matrix(MIDLOC_biomdata_lifeline$counts)
View(otu_table_lifeline)

#save rbiom object to a file
MIDLOCMicrobiome_lifeline_savepath <- 
write_biom(MIDLOC_biomdata_lifeline, MIDLOCMicrobiome_lifeline_savepath)

#adding metaphlan to the MD data and subset
# setting library paths MD data
MIDLOC_biompath_MD <- 
MIDLOC_metaphlan_tables_MDpath <- 

# import biom files into R
MIDLOC_biomdata_MD <- as_rbiom(MIDLOC_biompath_MD)
MIDLOC_biomdata_MD

#adding otu names towards the otu's
#import tsv files
MIDLOC_metaphlan_tables_MD <- data.table::fread(MIDLOC_metaphlan_tables_MDpath)
View(MIDLOC_metaphlan_tables_MD)

#change the otu_1 numbering towards the otu names
MIDLOC_biomdata_MD$otus <- MIDLOC_metaphlan_tables_MD$V1
otu_table_MD <- as.matrix(MIDLOC_biomdata_MD$counts)
View(otu_table_MD)

#subsetting to sample only used for the MIDLOC microbiome analysis
metadata <- MIDLOC_biomdata_MD$metadata
T0_samples <- metadata$.sample[grepl('T0$', metadata$.sample)]
MIDLOC_biomdata_MD_subset <- subset(MIDLOC_biomdata_MD, MIDLOC_biomdata_MD$metadata$.sample %in% T0_samples)
MIDLOC_biomdata_MD_subset

#save rbiom object to a file
MIDLOCMicrobiome_MD_savepath <- 
write_biom(MIDLOC_biomdata_MD_subset, MIDLOCMicrobiome_MD_savepath)
