##############################################
#    Ovine Lymph Node Liver Fluke RNA-seq     #
#     Mapping sheep to human orthologs        #
###############################################


# Author: Amalia Naranjo Lucena
# Last updated on: June 2020


############################################
# 1 Load and/or install required packages #
############################################

library(here)
library(devtools)
library(tidyverse)
library(magrittr)
library(tidyr)


############################################
# 2 Working directory, data, and time zone #
############################################

# Check working directory
here()

# Define variables for subdirectories
sheep <- here("sheep_genes.csv")
human <- here("human_genes.csv")
orthologs <- here("orthologs.csv")
Allgenes <- here("Tables/Infected_AllGenes.csv")
tablesDir <- here("Tables/")


# Set time zone
Sys.setenv(TZ = "Europe/London")

#######################
# 3 Importing files  #
#######################

#Import files

sheep_df <- read_csv(sheep)
human_df <- read_csv(human)
orthologs_df <- read_csv(orthologs)
Allgenes_df <- read.csv(Allgenes)

#Check files

sheep_df
human_df
orthologs_df
Allgenes_df


##############################
#  4 Cleaning sheep csv file #
##############################

# Sheep dataframe

sheep_df %<>%  
  dplyr::select(-c(gene_sheep, sheep_description)) %>%
  dplyr::rename(sheep_gene_symbol = gene_sheep_name) %>%
  dplyr::rename(Sheep_EntrezID = xref_sheep) %>%
  dplyr::filter(str_detect(Sheep_EntrezID, 'ncbi.nlm.nih.gov/gene/'))

# Check dataframe

View(sheep_df)

# Cleaning entrezid

sheep_df$Sheep_EntrezID %<>%
  str_remove("http://www.ncbi.nlm.nih.gov/gene/") %>%
  as.numeric()

# Check dataframe

View(sheep_df)

##############################
#  5 Cleaning human csv file #
##############################

# Human dataframe

human_df %<>%  
  dplyr::select(-c(gene_human, human_description)) %>%
  dplyr::rename(human_gene_symbol = gene_human_name) %>%
  dplyr::rename(Human_EntrezID = xref_human) %>%
  dplyr::filter(str_detect(Human_EntrezID, 'ncbi.nlm.nih.gov/gene/'))

# Check dataframe

View(human_df)

# Cleaning entrezid

human_df$Human_EntrezID %<>%
  str_remove("http://www.ncbi.nlm.nih.gov/gene/") %>%
  as.numeric()

# Check dataframe

View(human_df)


##################################
#  5 Cleaning orthologs csv file #
##################################


# Orthologs dataframe

orthologs_df %<>%  
  dplyr::select(-c(og, gene_sheep, gene_human, og_description)) %>%
  dplyr::rename(human_gene_symbol = gene_human_name) %>%
  dplyr::rename(sheep_gene_symbol = gene_sheep_name) 
  
# Check dataframe

View(orthologs_df)


########################
#  5 Cleaning Allgenes #
#######################

View(Allgenes_df)

# Rename columns

Allgenes_df %<>%
  dplyr::rename(sheep_gene_symbol = gene) %>%
  dplyr::rename(Sheep_EntrezID = EntrezID)
  

##########################
#  5 Matching sheep IDs  #
##########################


# Joining Allgenes to sheep.csv by BOTH EntrezID and gene symbol

Allgenes_df %>% 
  dplyr::left_join(sheep_df, by = c("Sheep_EntrezID", "sheep_gene_symbol")) %>%
  dplyr::select(Sheep_EntrezID, sheep_gene_symbol, everything()) -> sheep_join


# Check dimensions
dim(Allgenes_df)
dim(sheep_join)


# Get duplicated IDs

Duplicates <- c(sheep_join$Sheep_EntrezID[duplicated(sheep_join$Sheep_EntrezID)])
dplyr::filter(sheep_join, Sheep_EntrezID %in% Duplicates)

# Keep distinct IDs

sheep_join %<>%
  dplyr::distinct(.keep_all = FALSE)
 
View(sheep_join)

#Check dimensions

dim(Allgenes_df)
dim(sheep_join)

####################################
#  6 Join sheep_join to orthologs  #
####################################


# Joining sheep_join to to orthologs_df by sheep_gene_symbol

sheep_join %<>% 
  dplyr::left_join(orthologs_df, by = c("sheep_gene_symbol")) %>% 
  dplyr::select(Sheep_EntrezID, sheep_gene_symbol, human_gene_symbol, everything())


# Check dimensions
dim(Allgenes_df)
dim(sheep_join)


# Get duplicated IDs

Duplicates <- c(sheep_join$Sheep_EntrezID[duplicated(sheep_join$Sheep_EntrezID)])
View(dplyr::filter(sheep_join, Sheep_EntrezID %in% Duplicates))

# Keep distinct IDs and remove missing human gene symbols (NA)

sheep_join %<>%
  dplyr::distinct(.keep_all = FALSE) %>% 
  dplyr::filter(!is.na(human_gene_symbol))         
  

# Keep only rows that show EntrezID only once (trying to keep ONLY one to one orthologs)

Duplicates <- unique(sheep_join$Sheep_EntrezID[duplicated(sheep_join$Sheep_EntrezID)])

sheep_join %<>%
  dplyr::filter(!(Sheep_EntrezID %in% Duplicates)) 

# Check df

View(sheep_join)

# 7 Separate column with 2 IDs

sheep_join %<>% 
  separate(human_gene_symbol, c("Human_ID1", "Human_ID2"), sep = ";")

View(sheep_join)

# Save csv file

write_csv(sheep_join, file.path(paste0(tablesDir, "Orthologs_IPA.csv")),
          col_names = TRUE)


##########################
# 7 Save R session info #
##########################

devtools::session_info()

