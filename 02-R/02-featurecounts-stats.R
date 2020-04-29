###############################################
#    Ovine Lymph Node Liver Fluke RNA-seq     #
#       featurecounts statistics              #
###############################################

# Author: Amalia Naranjo Lucena
# Date: April 2020


############################################
# 01 Load and/or install required packages #
############################################

library(here)
library(tidyverse)
library(magrittr)
library(devtools)

# Uncomment functions below to install packages in case you don't have them
# CRAN packages
#install.packages("here")
#install.packages("tidyverse")
#install.packages("devtools")

######################################
# 02 Working directory and time zone #
######################################


# Check working directory
here()

# Define variables for subdirectories
summaryDir <- here("feature_Counts_summary")
tablesDir <- here("Tables/")

# Set time zone
Sys.setenv(TZ = "Europe/London")

##################################################
# 03 Import and tidy featureCounts summary files #
##################################################

# List files
summ_files <- list.files(summaryDir, full.names = TRUE)

# Import, tidy and export data:
map_dfc(summ_files, ~ read_tsv(.x)) %>%
  dplyr::rename(S = Status) %>%
  dplyr::select(-contains("Status")) %>%
  dplyr::rename_at(vars(contains("/")),
                   list(~str_remove(., "/home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/N\\d+_S\\d\\d/"))) %>%
  dplyr::rename_at(vars(contains("_")),
                   list(~str_remove(., "_Aligned.out.bam"))) %>%
  column_to_rownames(var = "S") %>%
  as.matrix() %>%
  t() %>%
  write.csv(file.path(paste0(tablesDir, "stats-featureCounts.csv")))


##########################
# 04 Save R session info #
##########################

devtools::session_info()



