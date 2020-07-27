###############################################
#    Ovine Lymph Node Liver Fluke RNA-seq     #
#       ngsShoRT filtering statistics         #
###############################################

# Author: Amalia Naranjo Lucena
# Date: July 2020

############################################
# 01 Load and/or install required packages #
############################################

library(here)
library(tidyverse)
library(magrittr)
library(devtools)
library(tidyr)

######################################
# 02 Working directory and time zone #
######################################


# Check working directory
here()
tablesDir <- here("Tables/")

# Define variables for subdirectories

DEgenes <- here("Infected_DEgenes.csv")
Davidunmapped <- here("DAVID_unmappedDEgenes.csv")


# Set time zone
Sys.setenv(TZ = "Europe/London")
 

####################
# 03 Import files #
####################

#Import files

DEgenes_df <- read_csv(DEgenes)
unmapped_df <- read_csv(Davidunmapped)

#Check files
DEgenes_df
unmapped_df

#######################
# 04 Join data frames #
#######################

#Join AllDEgenes and unmapped by EntrezID to remove unmapped genes

DEgenes_df %>% 
  dplyr::anti_join(unmapped_df, by = c("EntrezID"))  -> join
              
View(join)

# Save list of genes without unmapped

join %>%
  write_csv(file.path(paste0(tablesDir, "AllDAVIDgenes.csv")),
          col_names = TRUE)

############################################
# 05 filter values of logFC > 1 , or < -1 #
############################################

# LogFC Value were adjusted so that 3000 or less genes were kept
# No need to filter by FDR because DE genes was already filtered for
# FDR < 0.05 when it was created

join %>%
  dplyr::filter(abs(logFC) >= 0.49) -> david_df

# Save csv file

david_df %>%
  write_csv(file.path(paste0(tablesDir, "CutoffDAVIDgenes.csv")),
            col_names = TRUE)

##########################
# 06 Save R session info #
##########################

devtools::session_info()

