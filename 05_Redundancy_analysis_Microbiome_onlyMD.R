#Loading packages
library(rbiom)
library(phyloseq)
library(vegan)
library(dplyr)

# setting library paths
MIDLOCMicrobiom_biompath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/MD data/Diet/04_subsetted_data/04_subsetted_MIDLOCMicrobiome.biom"

# import biom file into R
MIDLOCMicrobiom_biomdata <- as_rbiom(MIDLOCMicrobiom_biompath)
MIDLOCMicrobiom_biomdata

#add new metadata
MIDLOC_metadatapath <- "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/metdata.tsv"
#add metadata to R
metadata<-readr::read_tsv(MIDLOC_metadatapath)
View(metadata)
metadata <- metadata[1:30,]

#could also rename if there are a lot of columns (sample is the name of the first column)
metadata <- rename(metadata, .sample = sample)

#renaming the column names
MIDLOCMicrobiom_biomdata$samples
MIDLOCMicrobiom_biomdata$samples[1:9]<- paste0("MIDLOC0", 10:18)
MIDLOCMicrobiom_biomdata$samples[10:16]<- paste0("MIDLOC00", 1:7)
MIDLOCMicrobiom_biomdata$samples[17]<- paste0("MIDLOC00", 9)
MIDLOCMicrobiom_biomdata$samples[18:28]<- paste0("MIDLOC0", 19:29)
MIDLOCMicrobiom_biomdata$samples[29:30]<- paste0("MIDLOC0", 31:32)
MIDLOCMicrobiom_biomdata$samples

#add metadata to the biom file
MIDLOCMicrobiom_biomdata$metadata <- metadata
MIDLOCMicrobiom_biomdata

#cleaning with janitor
clean_metadata <- janitor::clean_names(MIDLOCMicrobiom_biomdata$metadata)
colnames(clean_metadata) <- gsub("^pat_", "", colnames(clean_metadata))
clean_metadata <- rename(clean_metadata, .sample = sample)
View(clean_metadata)
clean_metadata$laxantia_21 <- NULL
names(clean_metadata)[10] <- "diabetes"
names(clean_metadata)[11] <- "diabetes_code"
names(clean_metadata)[14] <- "laxatives"

clean_metadata$laxatives_code <- ifelse(clean_metadata$laxatives == "yes", 1,0)

MIDLOCMicrobiom_biomdata$metadata <- clean_metadata
MIDLOCMicrobiom_biomdata

#create an otu_table 
otu_table <- as.matrix(MIDLOCMicrobiom_biomdata$counts)

#creating the dbRDA
MIDLOC_rda <- dbrda(t(otu_table) ~ MIDLOCMicrobiom_biomdata$metadata$age+MIDLOCMicrobiom_biomdata$metadata$sex_code +MIDLOCMicrobiom_biomdata$metadata$nmdas+MIDLOCMicrobiom_biomdata$metadata$diabetes_code +MIDLOCMicrobiom_biomdata$metadata$heteroplasmy+ MIDLOCMicrobiom_biomdata$metadata$bmi+MIDLOCMicrobiom_biomdata$metadata$bloating_score+ MIDLOCMicrobiom_biomdata$metadata$diarrhea_score+ MIDLOCMicrobiom_biomdata$metadata$constipation_score + MIDLOCMicrobiom_biomdata$metadata$laxatives_code + MIDLOCMicrobiom_biomdata$metadata$probiotica,
                    data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray")
MIDLOC_rda_nmdas <- dbrda(t(otu_table) ~ MIDLOCMicrobiom_biomdata$metadata$nmdas, 
                          data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray")

anova(MIDLOC_rda, permutations = 9999, by = "margin")
#checking multicollinearity
vif_values <- vegan::vif.cca(MIDLOC_rda)
vif_values
vif_data <- as.data.frame(vif_values)


#plot RDA
#plot the data in a figure
library(tidyverse)
library(broom)

#diagnostic plots
Residuals <- residuals(MIDLOC_rda)
Fitted_values <- fitted(MIDLOC_rda)


par(mfrow = c(2,2))
plot(Fitted_values, Residuals, main = "Residuals vs Fitted")
lines(lowess(Fitted_values, Residuals), col = "red", lwd = 2)
qqnorm(residuals(MIDLOC_rda))
qqline(Residuals, col = "red")

#checking the explained variance with PERMANOVA
dist_matrix <- vegdist(t(otu_table), method = "bray")

#looping all the variables seperated in a rda
variables_RDA <- colnames(MIDLOCMicrobiom_biomdata$metadata)
variables_RDA <- variables_RDA[-c(1,2,5,6,7,11,13,16,18,20,21)]

RDAresults <- data.frame(
  variable = character(),
  adj.r.squared = numeric(),
  r.squared = numeric(),
  inertia = numeric(),
  Pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (var in variables_RDA) {
  # Prepare formula for dbRDA, e.g., dbRDA(dist_matrix ~ variable)
  formula <- as.formula(paste("dist_matrix ~", var))
  
  # Run dbRDA with that variable
  model <- dbrda(formula, data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray")  # or other distance
  
  # Extract adjusted explained dissimilarity (R2 adj)
  R2adj <- RsquareAdj(model)$adj.r.squared
  
  # Extract explained dissimilarity (R2)
  R2 <- RsquareAdj(model)$r.squared
  
  #Extract the inertia
  
  intertia <- model$CCA$eig
  
  # Extract p-value from anova
  set.seed(1)
  a <- anova(model, permutations = 99999)
  pval <- a$`Pr(>F)`[1]
  
  # Add to results
  RDAresults <- rbind(RDAresults, data.frame(
    Variable = var,
    adj.r.squared = R2adj,
    r.squared = R2,
    inertia = intertia,
    Pvalue = pval,
    stringsAsFactors = FALSE
  ))
}

print(RDAresults)

RDAresults$ExplainedPercentage <- RDAresults$r.squared*100
RDAresults <- RDAresults %>%
  mutate(Variable = factor(Variable, levels = Variable[order(r.squared, decreasing = FALSE)]))
RDAresults$p.adj <- p.adjust(RDAresults$Pvalue, method = "bonferroni")

RDAresults$color <- ifelse(RDAresults$p.adj < 0.05, "brown1", "cornflowerblue")



#plotting it 
library(ggplot2)
RDA_plot <- ggplot(RDAresults, aes(x = ExplainedPercentage, y = Variable, fill = color))+
  geom_col(alpha =  0.8)+
  scale_fill_identity(name = "Legend",              
                      labels = c("Significant", "Not Significant"),   
                      breaks = c("brown1", "cornflowerblue"), 
                      guide = "legend")+
  labs(x = "Explained dissimilarity (%)", y = "Study variables")+
  theme_minimal()+
  theme(axis.title.x = element_text(size = 13, face = "bold"),
        axis.title.y = element_text(size = 13, face = "bold"))
RDA_plot

library(patchwork)
RDA_plot + plot_annotation(
  title = "Explained dissimilarity of the gut microbiome in MD patients",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)

#Saving the data as TSV file
library(readr)
write_tsv(RDAresults, "C:/Users/Z141231/OneDrive - Radboudumc/Rdata/Microbiome data/Functional/05_RDA/RDAresults.tsv")

# analysing the best model with ordistep
mod0 <- dbrda(dist_matrix ~ 1, data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray") # Model with intercept only
mod1 <- dbrda(dist_matrix ~ nmdas+diabetes_code +heteroplasmy+bloating_score+ diarrhea_score+ constipation_score +laxatives_code + probiotica + age + sex_code +bmi, data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray") # Model with all explanatory variables
## With scope present, the default direction is "both"
set.seed(123)
mod <- ordiR2step(mod0, scope = formula(mod1), permutations = 9999, trace = TRUE)
mod
## summary table of steps
mod$anova

plot_mod <- ordiplot(mod, type = "points")
plot_mod


##taking bmi and constipation_score into account, since they had low p-value (0,05-0,1)
mod_ordistep <- dbrda(dist_matrix ~ nmdas + bmi + constipation_score, data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray")
mod_ordistep
anova(mod_ordistep, by = "margin")

plot_ordistep <- ordiplot(mod_ordistep, scaling = 1, type = "text")
plot_ordistep

#checking the corrections
nmdas_cor <- dbrda(dist_matrix ~ nmdas + bmi + constipation_score + Condition(age+sex), data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray")
nmdas_cor
anova(nmdas_cor, permutations = 9999, by = "margin")

