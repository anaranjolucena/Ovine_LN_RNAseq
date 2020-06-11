#All samples were tested for clustering using a dendogram
#Run part A previous to this, and then run this

rawCounts$samples$group %<>%
  str_replace("N(5|15|24|49|58|71|76|79)_S\\d\\d$", "Control") %>%
  str_replace("N(12|13|17|19|23|42|51|55|66|67|73)_S\\d\\d$", "Infected") %>%
  factor(levels = c("Control", "Infected"))

LN_dgelist <- DGEList(counts       = rawCounts$counts,
                      group        = rawCounts$samples$group,
                      #genes = gene_annot,
                      lib.size     = NULL,
                      norm.factors = NULL,
                      remove.zeros = FALSE)


normalized.counts=cpm(LN_dgelist)
transposed=t(normalized.counts) # transposes the counts matrix
distance=dist(transposed) # calculates distance
clusters=hclust(distance) # does hierarchical clustering

# Open a pdf file
pdf(file.path(paste0(imgDir, "Clusterdendrogram_allsamples.pdf")))

# Create a plot
plot(clusters, main = "Cluster dendrogram showing all samples") # plots the clusters as a dendrogram

# Close the pdf file
dev.off() 

