#Loading packages
library("biomformat")
library(rbiom)

# setting library paths MD data
MIDLOC_biompath_Lifeline <- 
MIDLOC_biompath_MD <- 

# import biom files into R
MIDLOC_biomdata_MD <- as_rbiom(MIDLOC_biompath_MD)
MIDLOC_biomdata_MD

MIDLOC_biomdata_lifeline <- as_rbiom(MIDLOC_biompath_Lifeline)
MIDLOC_biomdata_lifeline

#check the otu_table 
otu_table_MD <- as.matrix(MIDLOC_biomdata_MD$counts)
colSums(otu_table_MD)
View(otu_table_MD)

otu_table_lifeline <- as.matrix(MIDLOC_biomdata_lifeline$counts)
colSums(otu_table_lifeline)
View(otu_table_lifeline)

#checking shared otu's
length(intersect(rownames(otu_table_MD), rownames(otu_table_lifeline)))
length(setdiff(rownames(otu_table_MD), rownames(otu_table_lifeline)))
length(setdiff(rownames(otu_table_lifeline), rownames(otu_table_MD)))

#checking the properties
biom_list <- list(MIDLOC_biomdata_lifeline, MIDLOC_biomdata_MD)

#remove metadata until I have metadata for lifelines
MIDLOC_biomdata_MD$metadata <- NULL
View(MIDLOC_biomdata_MD$metadata)

#trying to merge the biom files
MIDLOC_microbiomedata <- biom_merge(MIDLOC_biomdata_MD, MIDLOC_biomdata_lifeline)
MIDLOC_microbiomedata

otu_table <- as.matrix(MIDLOC_microbiomedata$counts)
colSums(otu_table)
View(otu_table)

#checking the merged biom
biom_list <- list(MIDLOC_biomdata_lifeline, MIDLOC_biomdata_MD, MIDLOC_microbiomedata)

#checking taxonomy
View(MIDLOC_microbiomedata$taxonomy)
View(MIDLOC_biomdata_MD$taxonomy)
View(MIDLOC_biomdata_lifeline$taxonomy)

different_values_MD <- length(setdiff(MIDLOC_biomdata_MD$taxonomy$.otu, MIDLOC_microbiomedata$taxonomy$.otu ))
different_values_MD

different_values_lifeline <- length(setdiff(MIDLOC_biomdata_lifeline$taxonomy$.otu, MIDLOC_microbiomedata$taxonomy$.otu))
different_values_lifeline

#checking if rowSums
not_zero <- rowSums(otu_table) != 0

#checking the alpha diversity
ad_MIDLOCMicrobiome <- adiv_table(MIDLOC_microbiomedata, adiv = ".all")
View(ad_MIDLOCMicrobiome)
ad_MIDLOC_lifeline <- adiv_table(MIDLOC_biomdata_lifeline, adiv = ".all")
View(ad_MIDLOC_lifeline)
ad_MIDLOC_MD <- adiv_table(MIDLOC_biomdata_MD, adiv = ".all")
View(ad_MIDLOC_MD)

#lifelines
shared_values <- intersect(ad_MIDLOCMicrobiome$.diversity, ad_MIDLOC_lifeline$.diversity)
shared_values

different_values <- length(setdiff(ad_MIDLOC_lifeline$.diversity, ad_MIDLOCMicrobiome$.diversity ))
different_values

#MD
shared_values <- intersect(ad_MIDLOCMicrobiome$.diversity, ad_MIDLOC_MD$.diversity)
shared_values

different_values <- length(setdiff(ad_MIDLOC_MD$.diversity, ad_MIDLOCMicrobiome$.diversity ))
different_values

#save rbiom object to a file
MIDLOCMicrobiomesavepath <- 
write_biom(MIDLOC_microbiomedata, MIDLOCMicrobiomesavepath)
