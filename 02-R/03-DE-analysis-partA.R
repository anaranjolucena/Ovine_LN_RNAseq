###############################################
#    Ovine Lymph Node Liver Fluke RNA-seq     #
#       featurecounts statistics              #
###############################################


# Author: Amalia Naranjo Lucena
# Last updated on: April 2020


############################################
# 01 Load and/or install required packages #
############################################

library(here)
library(tidyverse)
library(magrittr)
library(devtools)
library(extrafont)
library(edgeR)
library(biobroom)
library(Cairo)
library(rtracklayer)



# Uncomment functions below to install packages in case you don't have them
# Bioconductor packages
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("edgeR")
#BiocManager::install("biobroom")
#BiocManager::install("rtracklayer")

# CRAN packages
#install.packages("here")
#install.packages("tidyverse")
#install.packages("Cairo")
#install.packages("extrafont")
#install.packages("remotes")
#remotes::install_github("romainfrancois/nothing")

##############################################
# 02 Working directory, fonts, and time zone #
##############################################

# Check working directory
here()

# Define variables for subdirectories
countsDir <- here("feature_Counts")
imgDir <- here("Figures/")
tablesDir <- here("Tables/")
ovinegtf <- here("GCF_002742125.1_Oar_rambouillet_v1.0_genomic.gtf")

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device - font_import() only once
#font_import()
loadfonts()


#########################################
# 03 Import featureCounts sense counts  #
#########################################

# Create a vector of file names
files <- list.files(path        = countsDir,
                    pattern     = "^N",
                    all.files   = TRUE,
                    full.names  = FALSE,
                    recursive   = FALSE,
                    ignore.case = FALSE)

files
length(files)

# Create a dataframe with raw counts for all samples
rawCounts <- readDGE(path         = countsDir,
                     files        = files,
                     header       = TRUE,
                     comment.char = "#",
                     columns      = c(1, 7))
names(rawCounts)
head(rawCounts$samples)
head(rawCounts$counts)


#################################
# 04 Clean column and row names #
#################################

# Edit sample names
colnames(rawCounts$counts) %<>%
  str_replace("_counts", "")

rownames(rawCounts$samples) %<>%
  str_replace("_counts", "")


# Correct row names in counts
rownames(rawCounts$counts) %<>%
  str_replace(",miRBase:.*", "") %>%
  str_replace("GeneID:", "")

# Check data frames
head(rawCounts$counts)
head(rawCounts$samples)


#############################################
# 05 Extract  gene annotation from GTF file #
#############################################


#ANNOTATION RELEASE NAME:	NCBI Ovis aries Annotation Release 103
#ANNOTATION EVIDENCE FREEZE DATE:	25-January-2019
#ANNOTATION RELEASE DATE:	06-February-2019
#ANNOTATION REPORT:	https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Ovis_aries/103
#ASSEMBLY NAME:	Oar_rambouillet_v1.0
#ASSEMBLY ACCESSION:	GCF_002742125.1

#Importing gtf and transform to data frame
rtracklayer::import(ovinegtf,
                    format = "gtf") %>%
  as.data.frame() -> gtf

#Check data frame
View(gtf)

#Edit data frame to keep only: gene, db_xref and gbkey
gtf %>% 
  dplyr::select(db_xref, gene, gbkey) %>%
  dplyr::filter(gbkey == "Gene") -> gtf_annot

#Check data frame
head(gtf_annot)

# Cleaning entrezid
gtf_annot$db_xref %<>%
  str_replace("GeneID:", "") %>%
  as.numeric()

#Check data frame
head(gtf_annot)

#Getting EntrezID from rawCounts
gene_annot <- as.data.frame(rownames(rawCounts$counts)) 
colnames(gene_annot) <- c("EntrezID")

#Check dataframe
head(gene_annot)

#Check dimensions
dim(gene_annot)
dim(gtf_annot)
length(unique(gtf_annot$db_xref))


#Change EntrezID column to numeric
gene_annot$EntrezID <- as.numeric(as.character(gene_annot$EntrezID))

#Join gene annotations and remove gbkey column
gene_annot %<>% 
  left_join(gtf_annot, by = c("EntrezID" = "db_xref"))  %>%
  dplyr::select(-gbkey)

#Check dimensions
dim(gene_annot)
dim(rawCounts$counts)

# Get duplicated ID
gene_annot$EntrezID[duplicated(gene_annot$EntrezID)]
dplyr::filter(gene_annot, EntrezID == "443058")

# EntrezID 443058 is duplicated but gene name is the same in both (SPP1)
#Delete only one duplicate
gene_annot %<>%
  dplyr::distinct(.keep_all = FALSE)

#Check dimensions
dim(gene_annot)

#Check rows by ID
dplyr::filter(gene_annot, EntrezID == "443058")


#########################################
# 06 Add sample information for DGElist #
#########################################

# Treatment group (avoid using underscores)
rawCounts$samples$group <- rownames(rawCounts$samples)
rawCounts$samples$group %<>%
  str_replace("N(5|15|24|49|58|71|76|79)_S\\d\\d$", "Control") %>%
  str_replace("N(12|13|17|19|23|42|51|55|66|67|73)_S\\d\\d$", "Infected") %>%
  factor(levels = c("Control", "Infected"))

# Check data frame
rawCounts$samples
levels(rawCounts$samples$group)

#####################
# 07 Create DGElist #
#####################

# Use newly assigned variables to create DGElist
LN_dgelist <- DGEList(counts       = rawCounts$counts,
                          group        = rawCounts$samples$group,
                          genes = gene_annot,
                          lib.size     = NULL,
                          norm.factors = NULL,
                          remove.zeros = FALSE)
#Checks
names(LN_dgelist)
dim(LN_dgelist)
head(LN_dgelist$counts)
head(LN_dgelist$samples)
head(LN_dgelist$genes)
levels(LN_dgelist$samples$group)

#Check that rownames and EntreID column in genes dataframe are the same
identical(rownames(LN_dgelist$genes), as.character(LN_dgelist$genes$EntrezID))

# Include addtional experimental information into DGElist
# Subject number
LN_dgelist$samples$subject <- rownames(LN_dgelist$samples)
LN_dgelist$samples$subject %<>%
  str_replace("_S\\d\\d$", "") %>%
  factor()

# Check data frame
LN_dgelist$samples
levels(LN_dgelist$samples$subject)

# Export sample information
LN_dgelist$samples %>%
  rownames_to_column(var = "sample_name") %>%
  write_csv(file.path(paste0(tablesDir, "LN-sample-info.csv")),
            col_names = TRUE)

################################################
# 08 Density plot: raw gene counts per library #
################################################

# Tidy DGElist and plot data

LN_dgelist %>%
  tidy() %>%
  ggplot() +
  geom_density(aes(x     = log10(count + 1),
                   group = sample), size = 0.1) +
  theme_bw(base_size = 12, base_family = "Arial") +
  ggtitle("Density of raw gene counts per sample") +
  ylab("Density of raw gene counts per sample") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_raw

density_raw

# Export image
ggsave("LN-density_plot_raw_counts.png",
       plot      = density_raw,
       device    = "png",
       limitsize = FALSE,
       dpi       = 300,
       path      = imgDir)

###########################################
# 09 Remove zero and lowly expressed tags #
###########################################

# Filter non expressed tags (all genes that have zero counts in all samples)
LN_no_zeros <- LN_dgelist[rowSums(LN_dgelist$counts) > 0, ]
dim(LN_no_zeros$counts)
head(LN_no_zeros$counts)
dim(LN_dgelist$counts)


# Filter lowly expressed tags, retaining only tags with
# more than 1 count per million in 8 or more libraries
# (8 libraries correspond to the smallest group)
LN_filt <- LN_no_zeros[rowSums(cpm(LN_no_zeros) > 1) >= 8, ]
dim(LN_filt$counts)
head(LN_filt$counts)

# Ouptut filtered counts
LN_filt$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "EntrezID") %>%
  write_csv(file.path(paste0(tablesDir, "LN-sense_filt_counts.csv")),
            col_names = TRUE)


##############################
# 10 Recompute library sizes #
##############################

LN_filt$samples$lib.size <- colSums(LN_filt$counts)
head(LN_filt$samples)
head(LN_dgelist$samples)

###########################################################################
# 11 Calculate normalisation factors using Trimmed Mean of M-values (TMM) #
###########################################################################

# With edgeR, counts are not transformed in any way after
# calculating normalisation factors
LN_norm <- calcNormFactors(LN_filt, method = "TMM")
head(LN_norm$samples)

##########################
# 12 Save R session info #
##########################

devtools::session_info()

#######################
# 13 Save .RData file #
#######################

# Detach all loaded packages (except base R packages):
require(nothing, quietly = TRUE)

# Save all environment objects:
save.image(file = "ovineLN-DE-Analysis.RData")


