#Loading packages
library(rbiom)
library(phyloseq)
library(readr)

# setting library paths
MIDLOC_biompath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/Raw_data/metaphlan_with_taxonomy.biom"
MIDLOC_metadatapath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/Raw_data/metadata_MD_complete.tsv"
MIDLOC_genefamiliespath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/Raw_data/merged_genefamilies.tsv"
MIDLOC_pathabundancepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/Raw_data/merged_pathabundance.tsv"
MIDLOC_pathcoveragepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/Raw_data/merged_pathcoverage.tsv"
MIDLOC_metaphlan_tables <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/Raw_data/merged_metaphlan_tables.tsv"
biomtree <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/TestDataRadboudumc/practicetestdata/rooted_tree.newick"

# import biom file into R (Rbiom(2.0.13 is needed to load the file))
MIDLOC_biomdata <- as_rbiom(MIDLOC_biompath)
MIDLOC_biomdata

# create a simple matrix (It is also possible to load it via phylosec)
#mtx <- matrix(
  #data     = floor(runif(24) * 1000), 
  #nrow     = 6, 
  #dimnames = list(paste0("OTU", 1:6), paste0("Sample", 1:4)) )
#mtx

# convert matrix to rbiom -----------------------
#MIDLOC_biomdata<- as_rbiom(biom = mtx)
#MIDLOC_biomdata

# import biom file into R
#file <- phyloseq::import_biom(MIDLOC_biompath)
#file

#MIDLOC_biomdata <- as_rbiom(biom = file)
#MIDLOC_biomdata

#add metadata to R
metadata<-readr::read_tsv(MIDLOC_metadatapath)
print(metadata)

#rename metadata columnnames
#colnames(metadata)<- c('.sample','age','sex', 'sex_code', 'height', 'weight', 'BMI', 'NMDAS', "Patient_type", "Phenotype", "Diabetes", "heteroplasmy", "Timepoint", 'Timepoint_code')

#rename biom samples names for the metadata
#new_names <- c("MIDLOC06T24")
#MIDLOC06T24_biomdata$samples <- new_names

#create an otu_table to check the counts
otu_table <- as.matrix(MIDLOC_biomdata$counts)

#could also rename if there are a lot of columns (sample is the name of the first column)
metadata <- rename(metadata, .sample = sample)

#add metadata to the biom file
MIDLOC_biomdata$metadata <- metadata
MIDLOC_biomdata

#add genefamilies + path abundance + path coverage into R
genefamiliestable <- readr::read_tsv(MIDLOC_genefamiliespath)
View(genefamiliestable)
pathabundancetable <- readr::read_tsv(MIDLOC_pathabundancepath)
View(pathabundancetable)
pathcoveragetable <- readr::read_tsv(MIDLOC_pathcoveragepath)
View(pathcoveragetable)
metaphlan_tables <- readr::read_tsv(MIDLOC_metaphlan_tables)
View(metaphlan_tables)

#add tree to the biom file
biom$tree <- biomtree
biom

#save rbiom object to a file
MIDLOC_biomdata_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/01_Loaded_data/01_loaded_MIDLOC.biom"
write_biom(MIDLOC_biomdata, MIDLOC_biomdata_savepath)

#save tsv object to a file
genefamiliestable_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/01_Loaded_data/01_loaded_MIDLOCgenefamiliestable.tsv"
pathabundancetable_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/01_Loaded_data/01_loaded_MIDLOCpathabundancetable.tsv"
pathcoveragetable_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/01_Loaded_data/01_loaded_MIDLOCpathcoveragetable.tsv"
metaphlan_tables_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/01_Loaded_data/01_loaded_MIDLOCmetaphlan_tables.tsv"

write_tsv(genefamiliestable, genefamiliestable_savepath)
write_tsv(pathabundancetable, pathabundancetable_savepath)
write_tsv(pathcoveragetable, pathcoveragetable_savepath)
write_tsv(metaphlan_tables, metaphlan_tables_savepath)

#Loading multiqc report
multiqc_stats <- readr::read_tsv("C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/Raw_data/merged_multiqc_stats.tsv")
mean(multiqc_stats$`Input Reads`)
sd(multiqc_stats$`Input Reads`)
mean(multiqc_stats$`Clean Reads %`)
mean(multiqc_stats$`Clean Reads`)
sd(multiqc_stats$`Clean Reads`)
