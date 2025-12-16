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

#normalize by the root since the colsum is not comparable
library(slam)
feature_table <- MIDLOC_microbiomdata_pathabundance$counts
col_sums <- col_sums(feature_table)
features_table_norm <- feature_table
features_table_norm$v <- features_table_norm$v/col_sums[feature_table$j]
col_sums(features_table_norm)
MIDLOC_microbiomdata_pathabundance$counts <- features_table_norm

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

#adding metadata (can be placed here until better position)
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

#save rbiom object to a file TSS 1
MIDLOC_biomdata_savepath_TSS1 <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/04_Cleaned/04_cleaned_MIDLOC_TSS1.biom"
write_biom(MIDLOC_microbiomdata_pathabundance, MIDLOC_biomdata_savepath_TSS1)

#checking the DAA
Diff_pathway_table <- rbiom::taxa_stats(MIDLOC_microbiomdata_pathabundance, rank = 1, taxa = 0.0000001, stat.by = "group", test = "wilcox", unc = "drop")
View(Diff_pathway_table)

##From other script
##import tsv files
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
pathabundancetable_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/02_cleaned_data/02_cleaned_MIDLOCpathabundancetable.tsv"
pathcoveragetable_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/02_cleaned_data/02_cleaned_MIDLOCpathcoveragetable.tsv"

library(readr)
write_tsv(cleaned_pathabundancetable, pathabundancetable_savepath)
write_tsv(cleaned_pathcoveragetable, pathcoveragetable_savepath)


#adding protein name and go annotation to the gene familiestable
#example from the uniprot site
#import tsv files
genefamiliestable <- data.table::fread(MIDLOC_genefamiliespath)
View(genefamiliestable)
list_uniref90_codes <- cleaned_genefamiliestable$`# Gene Family`
list_uniref90_codes_1000 <- list_uniref90_codes[1:1000]
library(httr)

#Trying to loop the uniref90 codes to protein id via the uniprot site
isJobReady <- function(jobId) {
  pollingInterval = 5
  nTries = 20
  for (i in 1:nTries) {
    url <- paste0("https://rest.uniprot.org/idmapping/status/", jobId)
    r <- GET(url = url, accept_json())
    status <- content(r, as = "parsed")
    if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
      return(TRUE)
    }
    if (!is.null(status[["messages"]])) {
      print(status[["messages"]])
      return(FALSE)
    }
    Sys.sleep(pollingInterval)
  }
  return(FALSE)
}

getResultsURL <- function(redirectURL) {
  if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
    url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
  } else {
    url <- gsub("/results/", "/results/stream/", redirectURL)
  }
  return(url)
}

chunk_vector <- function(vec, chunk_size) {
  split(vec, ceiling(seq_along(vec) / chunk_size))
}

chunks <- chunk_vector(list_uniref90_codes, 1000)

all_results <- NULL

#since I hit the limit of uniprot so it needs multiple runs
start_chunk <- 1558

for (i in seq_along(chunks)) {
  if (i < start_chunk) {
    cat("Skipping chunk", i, "already processed.\n")
    next
  }
  chunk_ids <- chunks[[i]]
  cat("Processing chunk", i, "with", length(chunk_ids), "IDs...\n")
  
  files <- list(
    ids = paste(chunk_ids, collapse = ","),
    from = "UniRef90",
    to = "UniProtKB"
  )
  
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files, encode = "multipart", accept_json())
  submission <- content(r, as = "parsed")
  jobId <- submission[["jobId"]]
  
  if (is.null(jobId)) {
    warning("No job ID returned for chunk, skipping this chunk.")
    next
  }
  
  if (isJobReady(jobId)) {
    url <- paste0("https://rest.uniprot.org/idmapping/details/", jobId)
    r_details <- GET(url = url, accept_json())
    details <- content(r_details, as = "parsed")
    
    results_url <- getResultsURL(details[["redirectURL"]])
    results_url <- paste0(results_url, "?format=tsv")
    
    r_results <- GET(results_url)
    tsv_text <- content(r_results, as = "text")
    
    # Read the TSV results into a data frame
    chunk_results <- read_tsv(tsv_text, progress = FALSE)
    
    write_csv(chunk_results, paste0("C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/Gene_family_chunks/chunk_results_", i, ".csv"))
    rm(chunk_results)
    gc()
    
  } else {
    warning("Job did not complete for chunk with IDs: ", paste(head(chunk_ids, 5), collapse = ", "), " ...")
  }
  Sys.sleep(1)
}

# View combined results for all chunks
folder_path <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/Gene_family_chunks"
chunk_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
list_of_chunkfiles <- lapply(chunk_files, read_csv)
list_of_chunkfiles <- lapply(list_of_chunkfiles, function(df) {
  if ("Length" %in% colnames(df)) {
    # Only convert if it's not already numeric
    if (!is.numeric(df$Length)) {
      df$Length <- as.numeric(df$Length)
    }
  }
  return(df)
})
all_results <- bind_rows(list_of_chunkfiles)
View(all_results)

no_uncharacterizedproteins <- all_results[!grepl("Uncharacterized protein", all_results$`Protein names`, ignore.case = TRUE),]

#remove all the duplicate uniref90 and protein names, since the 
colnames(no_uncharacterizedproteins)[5] <- "Protein_names" 
cleaned_proteins <- no_uncharacterizedproteins %>%
  distinct(From, Protein_names, .keep_all = TRUE)

cleaned_proteins_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/02_cleaned_data/protein_names.tsv"
write_tsv(cleaned_proteins, cleaned_proteins_savepath)

colnames(cleaned_proteins)[1] <- "UniRef90_clean"
colnames(genefamiliestable)[1] <- "UniRef90"
only_proteinname <- cleaned_proteins[c("UniRef90_clean","Protein_names")]
genefamiliestable <- genefamiliestable %>%
  mutate(UniRef90_clean = sub("\\|.*", "", UniRef90))%>%
  select(UniRef90_clean, everything())
genefamiliestable_with_proteinnames <- genefamiliestable %>%
  inner_join(only_proteinname, by = "UniRef90_clean")%>%
  select(Protein_names, everything())
genefamiliestable_with_proteinnames$UniRef90_clean <- NULL 

#testing with a smaller data.frame
sample_genefamiliestable <- genefamiliestable[1:500,]
sample_genefamiliestable <- sample_genefamiliestable %>%
  mutate(UniRef90_clean = sub("\\|.*", "", UniRef90))%>%
  select(UniRef90_clean, everything())

cleaned_proteins <- readr::read_tsv("C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/02_cleaned_data/protein_names.tsv")
genefamiliestable_with_proteinnames <- sample_genefamiliestable %>%
  inner_join(only_proteinname, by = "UniRef90_clean")%>%
  select(Protein_names, everything())

#clean the column names so there is only the sample_number
colnames(genefamiliestable_with_proteinnames)[-c(1,2)] <- sub("_.*", "",colnames(genefamiliestable_with_proteinnames)[-c(1,2)])

#checking the data.frame
non_zero_rows <- genefamiliestable_with_proteinnames %>%
  filter(if_all(-c(1, 2), ~ . != 0))
filtered_rows <- genefamiliestable_with_proteinnames %>%
  filter(rowSums(select(., -c(1, 2)) == 0) <= 10)


#save tsv object to a file
genefamiliestable_savepath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/02_cleaned_data/02_cleaned_MIDLOCgenefamiliestable.tsv"

write_tsv(genefamiliestable_with_proteinnames, genefamiliestable_savepath)


# check the distribution
hist(otu_table, ylim = range(0, 1000))
log_otu_table <- log10(otu_table + 1e-5)
hist(log_otu_table, ylim = range(0, 25000))


hist_pathabundance <-column_to_rownames(cleaned_pathabundancetable, var = "# Pathway")
hist_pathabundance <- as.matrix(hist_pathabundance)
hist(hist_pathabundance, breaks = 45, ylim = range(0,25))

log_pathabundance <- log10(hist_pathabundance + 1e-5)
hist(log_pathabundance)
