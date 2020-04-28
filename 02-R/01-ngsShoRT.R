###############################################
#    Ovine Lymph Node Liver Fluke RNA-seq     #
#       ngsShoRT filtering statistics         #
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
ngsshort <- here("ngsshort.txt")
tablesDir <- here("Tables/")

# Set time zone
Sys.setenv(TZ = "Europe/London")


####################################
# 03 Import ngsShoRT summary files #
####################################

#Importing ngsshort file
ngstab <- read_table2(ngsshort)


#######################
# 04 Clean data frame #
#######################

#Retrieve names of columns
colnames(ngstab)

# Remove unnecessary columns
ngstab %<>%
  dplyr::select(-contains("X"))

# Remove extra symbols and convert data to numeric
ngstab$Percent_removed %<>% 
  str_replace("\\(", "") %>% 
  str_replace("\\)", "") %>% 
  str_replace("%", "") %>% 
  as.numeric()

#Check dataframe
ngstab


##################
# 05 Export data #
##################

# Export tidy summary stats
write_csv(ngstab,
          path = file.path(paste0(tablesDir, "ngsShort_stats.csv")),
          col_names = TRUE)

############################
# 06 R session information #
############################

devtools::session_info()
