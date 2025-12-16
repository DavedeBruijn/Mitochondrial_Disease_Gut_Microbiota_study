#Loading packages
library("biomformat")
library(rbiom)

# setting library paths MD data
MIDLOC_biompath_Lifeline <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/02_subsetted_andBiom/02_subsetted_andBiom_pathabundance_lifeline.biom"
MIDLOC_biompath_MD <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/02_subsetted_andBiom/02_subsetted_andBiom_pathabundance.biom"
MIDLOC_biompath_T1D <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/T1D data/02_metaphlan_subsetting/02_subsetted_andBiom_pathabundance_T1D.biom"

# import biom files into R
MIDLOC_biomdata_MD <- as_rbiom(MIDLOC_biompath_MD)
MIDLOC_biomdata_MD

MIDLOC_biomdata_lifeline <- as_rbiom(MIDLOC_biompath_Lifeline)
MIDLOC_biomdata_lifeline

MIDLOC_biomdata_T1D <- as_rbiom(MIDLOC_biompath_T1D)
MIDLOC_biomdata_T1D

#check the otu_table 
otu_table_MD <- as.matrix(MIDLOC_biomdata_MD$counts)
colSums(otu_table_MD)
View(otu_table_MD)
mean(colSums(otu_table_MD))

otu_table_lifeline <- as.matrix(MIDLOC_biomdata_lifeline$counts)
colSums(otu_table_lifeline)
View(otu_table_lifeline)
mean(colSums(otu_table_lifeline))

otu_table_T1D <- as.matrix(MIDLOC_biomdata_T1D$counts)
colSums(otu_table_T1D)
View(otu_table_T1D)
mean(colSums(otu_table_T1D))

#checking shared otu's
length(intersect(rownames(otu_table_MD), rownames(otu_table_lifeline)))
length(setdiff(rownames(otu_table_MD), rownames(otu_table_lifeline)))
length(setdiff(rownames(otu_table_lifeline), rownames(otu_table_MD)))

#checking the properties
biom_list <- list(MIDLOC_biomdata_MD, MIDLOC_biomdata_lifeline, MIDLOC_biomdata_T1D)

#remove metadata until I have metadata for lifelines
MIDLOC_biomdata_MD$metadata <- NULL
View(MIDLOC_biomdata_MD$metadata)

#trying to merge the biom files
MIDLOC_microbiomedata <- biom_merge(MIDLOC_biomdata_MD, MIDLOC_biomdata_lifeline, MIDLOC_biompath_T1D)
MIDLOC_microbiomedata

otu_table <- as.matrix(MIDLOC_microbiomedata$counts)
colSums(otu_table)
View(otu_table)

#checking the merged biom
biom_list <- list(MIDLOC_biomdata_lifeline, MIDLOC_biomdata_MD,MIDLOC_biomdata_T1D, MIDLOC_microbiomedata)

#checking taxonomy
View(MIDLOC_microbiomedata$taxonomy)
View(MIDLOC_biomdata_MD$taxonomy)
View(MIDLOC_biomdata_lifeline$taxonomy)
View(MIDLOC_biomdata_T1D$taxonomy)

different_values_MD <- length(setdiff(MIDLOC_biomdata_MD$taxonomy$.otu, MIDLOC_microbiomedata$taxonomy$.otu ))
different_values_MD

different_values_lifeline <- length(setdiff(MIDLOC_biomdata_lifeline$taxonomy$.otu, MIDLOC_microbiomedata$taxonomy$.otu))
different_values_lifeline

different_values_T1D <- length(setdiff(MIDLOC_biomdata_T1D$taxonomy$.otu, MIDLOC_microbiomedata$taxonomy$.otu))
different_values_T1D
#checking if rowSums
not_zero <- rowSums(otu_table) != 0
all(not_zero)

#checking the alpha diversity
ad_MIDLOCMicrobiome <- adiv_table(MIDLOC_microbiomedata, adiv = ".all")
View(ad_MIDLOCMicrobiome)
ad_MIDLOC_lifeline <- adiv_table(MIDLOC_biomdata_lifeline, adiv = ".all")
View(ad_MIDLOC_lifeline)
ad_MIDLOC_MD <- adiv_table(MIDLOC_biomdata_MD, adiv = ".all")
View(ad_MIDLOC_MD)
ad_MIDLOC_T1D <- adiv_table(MIDLOC_biomdata_T1D, adiv = ".all")
View(ad_MIDLOC_T1D)

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

#T1D
shared_values <- intersect(ad_MIDLOCMicrobiome$.diversity, ad_MIDLOC_T1D$.diversity)
shared_values

different_values <- length(setdiff(ad_MIDLOC_MD$.diversity, ad_MIDLOC_T1D$.diversity ))
different_values

#save rbiom object to a file
MIDLOCMicrobiomesavepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/03_Merged/03_Merged.biom"
write_biom(MIDLOC_microbiomedata, MIDLOCMicrobiomesavepath)
