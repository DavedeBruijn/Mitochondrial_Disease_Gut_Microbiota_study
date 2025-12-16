#Loading packages
library(phyloseq)
library(readr)
library(tibble)
library(dplyr)

#setting librarypaths
MD_pathcoverage_path <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/01_loaded_data/01_loaded_MIDLOCpathcoveragetable.tsv"
Control_pathcoverage_path <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Lifelines data/01_loaded/01_loaded_MIDLOCpathcoveragetable_lifeline.tsv"
T1D_pathcoverage_path <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/T1D data/01_loaded/01_loaded_MIDLOCpathcoveragetable.tsv"

#loading pathcoverage
pathcoveragetable_MD <- readr::read_tsv(MD_pathcoverage_path)
View(pathcoveragetable_MD)

pathcoveragetable_Control <- readr::read_tsv(Control_pathcoverage_path)
View(pathcoveragetable_Control)

pathcoveragetable_T1D <- readr::read_tsv(T1D_pathcoverage_path)
View(pathcoveragetable_T1D)

#not zero checking if rowSums
pathcoverage_rownames <- column_to_rownames(pathcoveragetable_MD, var = "# Pathway")
not_zero <- rowSums(pathcoverage_rownames) != 0
not_zero <- as.data.frame(not_zero)
all(not_zero$not_zero)
which(!not_zero$not_zero)

#combining columns
pathcoverage_total <- dplyr::left_join(pathcoveragetable_MD, pathcoveragetable_Control, pathcoveragetable_T1D, by = "# Pathway", keep = TRUE)

#filter the first row so we only have pathways (no pathways + genies/species) needed for diversity
filtered_pathcoveragetable_MD <- pathcoveragetable_MD[!grepl("\\|", pathcoveragetable_MD$`# Pathway`),]
filtered_pathcoveragetable_Control <- pathcoveragetable_Control[!grepl("\\|", pathcoveragetable_Control$`# Pathway`),]
filtered_pathcoveragetable_T1D <- pathcoveragetable_T1D[!grepl("\\|", pathcoveragetable_T1D$`# Pathway`),]

#not zero checking if rowSums
pathcoverage_rownames <- column_to_rownames(filtered_pathcoveragetable_Control, var = "# Pathway")
not_zero <- rowSums(pathcoverage_rownames) != 0
not_zero <- as.data.frame(not_zero)
all(not_zero$not_zero)
which(!not_zero$not_zero)

#calculating the pathcoverage per pathway.
pathcoverage_mean <- data.frame(pathway = pathcoveragetable_MD[,1], Mean = rowMeans(pathcoveragetable_MD[,-1]))
pathcoverage_mean_filtered <- data.frame(pathway = filtered_pathcoveragetable_MD[,1], Mean = rowMeans(filtered_pathcoveragetable_MD[,-1]))
pathcoverage_50_filtered <- pathcoverage_mean_filtered[pathcoverage_mean_filtered$Mean > 0.5, ]

pathcoverage_70 <- pathcoverage_mean[pathcoverage_mean$Mean > 0.7, ]
pathcoverage_50 <- pathcoverage_mean[pathcoverage_mean$Mean > 0.5, ]
pathcoverage_30 <- pathcoverage_mean[pathcoverage_mean$Mean > 0.3, ]


pathcoverage_mean_control <- data.frame(pathway = pathcoveragetable_Control[,1], Mean = rowMeans(pathcoveragetable_Control[,-1]))

pathcoverage_70_control <- pathcoverage_mean[pathcoverage_mean$Mean > 0.7, ]
pathcoverage_50_control <- pathcoverage_mean[pathcoverage_mean$Mean > 0.5, ]
pathcoverage_30_control <- pathcoverage_mean[pathcoverage_mean$Mean > 0.3, ]

pathcoverage_mean_T1D <- data.frame(pathway = pathcoveragetable_T1D[,1], Mean = rowMeans(pathcoveragetable_T1D[,-1]))

pathcoverage_70_T1D <- pathcoverage_mean[pathcoverage_mean$Mean > 0.7, ]
pathcoverage_50_T1D <- pathcoverage_mean[pathcoverage_mean$Mean > 0.5, ]
pathcoverage_30_T1D <- pathcoverage_mean[pathcoverage_mean$Mean > 0.3, ]

shared_pathways <- intersect(pathcoverage_30$X..Pathway, pathcoverage_30_control$X..Pathway)
shared_pathways

shared_pathways <- intersect(pathcoverage_30$X..Pathway, pathcoverage_30_T1D$X..Pathway)
shared_pathways

#getting the mean from a specific pathway
pathway <- pathcoveragetable_Control[pathcoveragetable_Control$`# Pathway`== "PENTOSE-P-PWY: pentose phosphate pathway", ]
rowMeans(pathway[,-1])

#save pathways with a coverage of 30 or higher
pathcoverage_30_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/00_SanityCheck/pathcoverage_30.tsv"
write_tsv(pathcoverage_30, pathcoverage_30_savepath)

#save only pathways with a coverage of 50 or higher
pathcoverage_50_filtered_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/00_SanityCheck/pathcoverage_50_filtered.tsv"
write_tsv(pathcoverage_50_filtered, pathcoverage_50_savepath)

pathcoverage_50_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/00_SanityCheck/pathcoverage_50.tsv"
write_tsv(pathcoverage_50, pathcoverage_50_savepath)
