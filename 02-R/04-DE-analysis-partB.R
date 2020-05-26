###############################################
#    Ovine Lymph Node Liver Fluke RNA-seq     #
#       featurecounts statistics              #
###############################################


# Author: Amalia Naranjo Lucena
# Last updated on: May 2020

############################################
# 14 Load and/or install required packages #
############################################

library(here)
library(edgeR)
library(devtools)
library(tidyverse)
library(magrittr)
library(biobroom)
library(ggridges)
library(Cairo)
library(extrafont)
library(ggrepel)
library(ggfortify)
library(statmod)


# Uncomment functions below to install packages in case you don't have them
#install.packages("ggridges")
#install.packages("ggrepel")
#install.packages("ggfortify")
#install.packages("statmod")


####################################################
# 15 Working directory, data, fonts, and time zone #
####################################################

# Check working directory
here()

# Load previously saved data
load("ovineLN-DE-Analysis.RData")

# Check variables for subdirectories
imgDir
tablesDir

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device
loadfonts()

##########################################################
# 16 Tidy DGElist for exploratory data analysis plotting #
##########################################################

tidy_LNNorm <- tidy(LN_norm, addSamples = TRUE)

# Turn sample info into factors for plotting
tidy_LNNorm$subject %<>%
  factor(levels = c("N5", "N12", "N13", "N15", "N17", 
                    "N19", "N23", "N24", "N42", "N49", 
                    "N51", "N55", "N67","N71", "N73", 
                    "N76", "N79"))

# Check factors
levels(tidy_LNNorm$group)
levels(tidy_LNNorm$subject)

# Check data frame
tidy_LNNorm



########################################################
# 17 Plot: density of filtered gene counts per library #
########################################################

ggplot(tidy_LNNorm, aes(x = log10(count + 1),
                            y = subject)) +
  scale_y_discrete(limits = rev(levels(tidy_LNNorm$subject))) +
  geom_density_ridges(aes(fill = group), alpha = 0.5) +
  scale_fill_manual("Condition", values = c("#b2b2b2", "#e06377")) +
  theme_bw(base_size = 12, base_family = "Arial") +
  ggtitle("Density of filtered gene counts per subject") +
  facet_grid(. ~ group) +
  ylab("Animal ID") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_norm

density_norm

# Export high quality image
ggsave("LN-density-filt.pdf",
       plot      = density_norm,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 9,
       units     = "in")


#########################################
# 18 Plot: principal component analysis #
#########################################

# Plot PCA using log2CPM filtered counts
pca <- prcomp(t(cpm(LN_norm, log = TRUE)),
              scale. = TRUE)
names(pca)

# Get PCA summary
summary_pca <- summary(pca)
summary_pca
names(summary_pca)

# Export PCA summary
summary_pca$importance %>%
  write.csv(file.path(paste0(tablesDir, "summary_pca.csv")))

# Plot component variance
data.frame(summary_pca$importance[2, ],
           stringsAsFactors = TRUE) %>%
  dplyr::rename("Variance" = "summary_pca.importance.2...") %>%
  rownames_to_column(var = "PC") %>%
  ggplot() +
  geom_col(aes(x = fct_inorder(PC), y = Variance)) +
  ylab("Proportion of Variance") +
  xlab("Principal Component") +
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) -> variance

variance

# Export high quality image
ggsave(filename  = "PCA_variance.pdf",
       plot      = variance,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 5,
       width     = 6,
       units     = "in")


# Edit sample names for plotting
edited_samples <- LN_norm$samples
rownames(edited_samples) %<>%
  str_remove("_S\\d\\d$")

edited_samples

# Plot PC1 versus PC2
pc1vspc2 <- autoplot(pca,
                     data = edited_samples,
                     colour = "group",
                     label = TRUE,
                     label.show.legend = FALSE,
                     label.repel = TRUE) +
  theme_bw(base_size = 12, base_family = "Arial")

pc1vspc2

# Export high quality image
ggsave(filename  = "pc1vspc2.pdf",
       plot      = pc1vspc2,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 5,
       width     = 7,
       units     = "in")

######################
# 19 Plot dendrogram #
######################

# Transpose the counts matrix, calculation of distance and hierarchical clustering

LN_dgelist %>%
  cpm() %>%
  t() %>%
  dist() %>%
  hclust() -> clusters

clusters

# Check what clusters is composed of
names(clusters)

# Check labels and remove _Sdd

clusters$labels

clusters$labels %<>%
  str_remove("_S\\d\\d")

clusters$labels

# Open a pdf file
pdf(file.path(paste0(imgDir, "Clusterdendrogram_samplesremoved.pdf")))

# Create a plot
plot(clusters, main = "Cluster dendrogram without N58 and N66") # plots the clusters as a dendrogram

# Close the pdf file
dev.off() 


###############################
# 20 Create the design matrix #
###############################

# Create a design matrix
design <- model.matrix(~group,
                       data = LN_norm$samples)

# Check design matrix
design

# Rename design matrix columns for simplicity
colnames(design) %<>%
  str_replace("group", "")

# Check design matrix
design
colnames(design)

#########################################
# 21 Estimate the dispersion parameters #
#########################################

# Common and trended dispersions are estimated with the
# Cox-Reid method and tagwise dispersions with the
# empirical Bayes method
LN_disp <- estimateDisp.DGEList(y       = LN_norm,
                                  design  = design,
                                  robust  = TRUE,
                                  verbose = TRUE)

names(LN_disp)

# Check the calculated dispersion
LN_disp$common.dispersion

# Check the calculated dispersion's square root,
# which corresponds to the biological coefficient of variation (BCV)
sqrt(LN_disp$common.dispersion)


################################
# 22 Plot: BCV and dispersions #
################################

# Create a dataframe with the dispersion values
names(LN_disp)

Disp <- as.data.frame(cbind(LN_disp$genes,
                            LN_disp$tagwise.dispersion,
                            LN_disp$common.dispersion,
                            LN_disp$trended.dispersion,
                            LN_disp$AveLogCPM))

colnames(Disp) %<>%
  str_replace("LN_disp\\$", "")

Disp %<>%
  dplyr::mutate(type_point = "Tagwise dispersion") %>%
  dplyr::mutate(type_hline = "Common dispersion") %>%
  dplyr::mutate(type_smooth = "Trended dispersion")

# Plot all dispersions
LN_BCV <- ggplot(Disp) +
  geom_point(aes(x = AveLogCPM,
                 y = sqrt(tagwise.dispersion),
                 fill = type_point),
             alpha = 0.5) +
  geom_hline(aes(yintercept = sqrt(common.dispersion),
                 colour = type_hline)) +
  geom_smooth(aes(x = AveLogCPM,
                  y = sqrt(trended.dispersion),
                  colour = type_smooth),
              linetype = 2) +
  scale_fill_manual("", values = c("black")) +
  scale_colour_manual("", values = c("red", "blue")) +
  theme_bw(base_size = 14, base_family = "Arial") +
  ggtitle("Estimated dispersions (NB model)") +
  xlab(expression(paste("Average ", log[2],"CPM"))) +
  ylab("Biological Coefficient of Variation")

LN_BCV

# Output high resolution plot
ggsave("LN_BCV.pdf",
       plot = LN_BCV,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 14,
       units     = "in")



##########################
# 23 Save R session info #
##########################

devtools::session_info()

#######################
# 24 Save .RData file #
#######################

# Detach all loaded packages (except base R packages):
require(nothing, quietly = TRUE)

# Save all environment objects:
save.image(file = "ovineLN-DE-Analysis.RData")


