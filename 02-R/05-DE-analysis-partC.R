###############################################
#    Ovine Lymph Node Liver Fluke RNA-seq     #
#       featurecounts statistics              #
###############################################


# Author: Amalia Naranjo Lucena
# Last updated on: May 2020

############################################
# 25 Load and/or install required packages #
############################################

library(here)
library(edgeR)
library(devtools)
library(tidyverse)
library(magrittr)
library(biobroom)
library(Cairo)
library(extrafont)
library(ggrepel)
library(ggfortify)
library(statmod)


# Uncomment functions below to install packages in case you don't have them
#install.packages("ggrepel")
#install.packages("ggfortify")
#install.packages("statmod")


####################################################
# 26 Working directory, data, fonts, and time zone #
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

####################
# 27 Fit GLM model #
####################

# Fit a quasi-likelihood negative binomial generalized log-linear model
# for each tag using the design matrix and calculated dispersion
LN_fitQL <- glmQLFit(y = LN_disp,
                       design = design,
                       robust = TRUE)

names(LN_fitQL)
colnames(LN_fitQL$design)


###########################
# 28 Plot: QL dispersions #
###########################

# Create a dataframe with the dispersion values
names(LN_fitQL)

DispQL <- as.data.frame(cbind(AveLogCPM = LN_fitQL$AveLogCPM,
                              deviance = LN_fitQL$deviance,
                              df.residual.zeros = LN_fitQL$df.residual.zeros,
                              var.prior = LN_fitQL$var.prior,
                              var.post = LN_fitQL$var.post))

DispQL %<>%
  dplyr::mutate(type_point = "Raw dispersion (NB)") %>%
  dplyr::mutate(type_point2 = "Squeezed EB dispersion") %>%
  dplyr::mutate(type_smooth = "Trended EB dispersion")

head(DispQL)


# Plot all dispersions
LN_BCVQL <- ggplot(DispQL) +
  geom_point(aes(x = AveLogCPM,
                 y = sqrt(sqrt(deviance/df.residual.zeros)),
                 fill = type_point),
             alpha = 0.5) +
  geom_point(aes(x = AveLogCPM,
                 y = sqrt(sqrt(var.post)),
                 colour = type_point2),
             alpha = 0.5) +
  geom_smooth(aes(x = AveLogCPM,
                  y = sqrt(sqrt(var.prior)),
                  colour = type_smooth),
              linetype = 2) +
  scale_fill_manual("", values = c("black")) +
  scale_colour_manual("", values = c("indianred4", "blue")) +
  theme_bw(base_size = 14, base_family = "Arial") +
  ggtitle("Estimated QL dispersions") +
  xlab(expression(paste("Average ", log[2],"CPM"))) +
  ylab("Quarter-Root Mean Deviance")

LN_BCVQL

# Output high resolution plot
ggsave("LN_BCVQL.pdf",
       plot = LN_BCVQL,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 14,
       units     = "in")

#####################################################
# 29 Determine differential expression by applying a #
#  Quasi-likelihood Tests                            #
#####################################################

# Test for differential expression,
# using the coefficient from LNfit$design

Infected.QL <- glmQLFTest(LN_fitQL, coef = "Infected")
testDE.Infected <- topTags(object        = Infected.QL,
                      n             = "inf",
                      adjust.method = "BH")

head(testDE.Infected$table)
dim(testDE.Infected$table)

### Output all genes tested for DE

testDE.Infected$table %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "true", "false")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, "Infected_AllGenes.csv")),
            col_names = TRUE)


### Filter genes considered DE (FDR < 0.05)

testDE.Infected$table %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(desc(logFC)) %>%
  as_tibble() -> InfectedDE

# Check dataframe

InfectedDE
tail(InfectedDE)

# Output
InfectedDE %>%
write_csv(file.path(paste0(tablesDir, "Infected_DEgenes.csv")),
          col_names = TRUE)

##################################
# 30 Plot: Volcano of DE genes  #
#################################


testDE.Infected$table %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "true", "false")) %>%
  ggplot(aes(x = logFC,
             y = -log10(FDR),
             colour = is.DE)) +
  geom_point(alpha = 0.4, size  = 1.75) +
  geom_label_repel(aes(label = if_else(abs(logFC) > 2.5 & FDR < 0.05,
                                       gene, NULL),
                       fontface = "italic",
                       family = "Arial",
                       size = 4,
                       fill = "gray",
                       alpha = 0.5),
                   show.legend = FALSE) +
  xlim(c(-5.5, 5.7)) +
  ylim(c(0, 7.5)) +
  theme_bw(base_size = 14, base_family = "Arial") +
  ggtitle("Infected DE genes") +
  xlab(expression(paste(log[2], " fold-change"))) +
  ylab(expression(paste(-log[10], " FDR"))) +
  scale_colour_manual("FDR \n< 0.05",
                      labels = c("False", "True"),
                      values = c('#41b6c4', '#225ea8')) -> Infected_Volcano

Infected_Volcano

# Output high resolution plot
ggsave("Infected_Volcano.pdf",
       plot = Infected_Volcano,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 14,
       units     = "in")


##########################
# 31 Save R session info #
##########################

devtools::session_info()

#######################
# 32 Save .RData file #
#######################

# Detach all loaded packages (except base R packages):
require(nothing, quietly = TRUE)

# Save all environment objects:
save.image(file = "ovineLN-DE-Analysis.RData")

