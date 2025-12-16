#Loading packages
library(rbiom)
library(phyloseq)
library(vegan)
library(dplyr)

# setting library paths
MIDLOCMicrobiom_biompath <- 

# import biom file into R
MIDLOCMicrobiom_biomdata <- as_rbiom(MIDLOCMicrobiom_biompath)
MIDLOCMicrobiom_biomdata

#renaming the rank
colnames(MIDLOCMicrobiom_biomdata$taxonomy) <- c(".otu", "pathway", "genus", "species")
MIDLOCMicrobiom_biomdata

#create an otu_table 
otu_table <- as.matrix(MIDLOCMicrobiom_biomdata$counts)
colSums(otu_table)

#normalize by the root since the colsum is not comparable
library(slam)
feature_table <- MIDLOCMicrobiom_biomdata$counts
col_sums <- col_sums(feature_table)
features_table_norm <- feature_table
features_table_norm$v <- features_table_norm$v/col_sums[feature_table$j]
col_sums(features_table_norm)
MIDLOCMicrobiom_biomdata$counts <- features_table_norm

#check the otu_table 
otu_table <- as.matrix(MIDLOCMicrobiom_biomdata$counts)
colSums(otu_table)
View(otu_table)

colnames(MIDLOCMicrobiom_biomdata$metadata)[14] <- "laxatives"

MIDLOCMicrobiom_biomdata$metadata$laxatives_score <- ifelse(MIDLOCMicrobiom_biomdata$metadata$laxatives == "yes", 1, 0)

#subsetting for only patients with MD
MIDLOCMicrobiom_biomdata <- subset(MIDLOCMicrobiom_biomdata, MIDLOCMicrobiom_biomdata$metadata$group %in% "MD")
MIDLOCMicrobiom_biomdata
View(MIDLOCMicrobiom_biomdata$metadata)

#create new otu_table
otu_table <- as.matrix(MIDLOCMicrobiom_biomdata$counts)
View(otu_table)

#creating the dbRDA
MIDLOC_rda <- dbrda(t(otu_table) ~ MIDLOCMicrobiom_biomdata$metadata$age+MIDLOCMicrobiom_biomdata$metadata$sex_code +MIDLOCMicrobiom_biomdata$metadata$nmdas+MIDLOCMicrobiom_biomdata$metadata$diabetes_code +MIDLOCMicrobiom_biomdata$metadata$heteroplasmy+ MIDLOCMicrobiom_biomdata$metadata$bmi+MIDLOCMicrobiom_biomdata$metadata$bloating_score+ MIDLOCMicrobiom_biomdata$metadata$diarrhea_score+ MIDLOCMicrobiom_biomdata$metadata$constipation_score + MIDLOCMicrobiom_biomdata$metadata$laxatives_score + MIDLOCMicrobiom_biomdata$metadata$probiotica,
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

#checking the explained variance 
dist_matrix <- vegdist(t(otu_table), method = "bray")

#looping all the variables seperated in a rda
variables_RDA <- colnames(MIDLOCMicrobiom_biomdata$metadata)
variables_RDA <- variables_RDA[-c(1,2,5,6,7,11,13,16,18,20)]

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
write_tsv(RDAresults, )

# analysing the best model with ordistep
mod0 <- dbrda(dist_matrix ~ 1, data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray") # Model with intercept only
mod1 <- dbrda(dist_matrix ~ nmdas+diabetes_code +heteroplasmy+bloating_score+ diarrhea_score+ constipation_score +laxatives + probiotica + age + sex_code +bmi, data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray") # Model with all explanatory variables
## With scope present, the default direction is "both"
set.seed(123)
mod <- ordiR2step(mod0, scope = formula(mod1), permutations = 9999, trace = TRUE)
mod
## summary table of steps
mod$anova

plot_mod <- ordiplot(mod, type = "points")
plot_mod


##taking bmi and constipation_score into account, since they had low p-value (0,05-0,1)
mod_ordistep <- dbrda(dist_matrix ~ nmdas + bmi + laxatives, data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray")
mod_ordistep
anova(mod_ordistep, by = "margin")

plot_ordistep <- ordiplot(mod_ordistep, scaling = 1, type = "text")
plot_ordistep

#checking the corrections
nmdas_cor <- dbrda(dist_matrix ~ nmdas + bmi + constipation_score + Condition(age+sex), data = MIDLOCMicrobiom_biomdata$metadata, distance = "bray")
nmdas_cor
anova(nmdas_cor, permutations = 9999, by = "margin")

