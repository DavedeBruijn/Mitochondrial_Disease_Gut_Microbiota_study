#Loading packages
library(rbiom)
library(phyloseq)
library(readr)

# setting library paths
MIDLOC_biompath <- 
MIDLOC_metadatapath <- 
MIDLOC_genefamiliespath <- 
MIDLOC_pathabundancepath <- 
MIDLOC_pathcoveragepath <- 
MIDLOC_metaphlan_tables <- 


# import biom file into R (Rbiom(2.0.13 is needed to load the file))
MIDLOC_biomdata <- as_rbiom(MIDLOC_biompath)
MIDLOC_biomdata

#add metadata to R
metadata<-readr::read_tsv(MIDLOC_metadatapath)
print(metadata)

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


#save rbiom object to a file
MIDLOC_biomdata_savepath <- 
write_biom(MIDLOC_biomdata, MIDLOC_biomdata_savepath)

#save tsv object to a file
genefamiliestable_savepath <- 
pathabundancetable_savepath <- 
pathcoveragetable_savepath <- 
metaphlan_tables_savepath <- 

write_tsv(genefamiliestable, genefamiliestable_savepath)
write_tsv(pathabundancetable, pathabundancetable_savepath)
write_tsv(pathcoveragetable, pathcoveragetable_savepath)
write_tsv(metaphlan_tables, metaphlan_tables_savepath)

#Loading multiqc report
multiqc_stats <- readr::read_tsv("filepath/merged_multiqc_stats.tsv")
mean(multiqc_stats$`Input Reads`)
sd(multiqc_stats$`Input Reads`)
mean(multiqc_stats$`Clean Reads %`)
mean(multiqc_stats$`Clean Reads`)
sd(multiqc_stats$`Clean Reads`)
