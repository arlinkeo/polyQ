# Main script

# Obtain data information of polyQ genes and brain regions
source("polyQ_region.R")
# Information is stored in "resources/polyQ.RData" and called whenever needed

# baseScript.R is called in every script to have basic function an variables needed

# Read microarray expression data and save as R-object 
source("brainExpr.R")
# Data saved as "resources/BrainExpr.RData"

# Select region-specific samples
# Save as binary vectors indicating for each sample whether it is within the specified region (1) or not (0)
source("sampleIDs.R")
# Sample selection saved as "resources/sampleIDs.RData"

# Plot pie charts of samples in HD region and samples in anatomical regions
source("pieChart_samples.R")

# Check number of samples of atomic structures that fall within the HD associated region (Coppen2016)
source("overlap_samples.R")

# Mean and variance of polyQ gene expression across regional-specific samples and all donors
source("avgExpr.R")
# Mean expr. is saved as "resources/avgExpr.RData"
# Mean and variance table is saved as heatmap in "expr_meanvarmedian.pdf"
# Color is assigned to the mean expression and saved as "resources/avgExprColor.RData"
# Legend bar of the expression is saved in "avgExpr.pdf"

# Differential expression of HD-associated region
# Analysis is done for each brain separately, then genes are selected if significant in at least 5 brains
# For the gene selection, functional enrichment analysis is doen for up- and down-regulated genes
source("diffExpr_HDregion.R")
# Differentially expressed genes are saved as "resources/diffgenesPerBrain.RData"
# A pdf is saved with  boxplots of the expression of each polyQ genes in- and outside the HD-region.

# Co-expression network of HD-associated region
# This script was ran on a server, because of memory requirements
source("coexpr_SCN_Coppen2016.R")
# correlation matrix saved as "meanCor_[region].RData"

# Plot histograms of correlation distribution of each polyQ gene and brain region
source("coexpr_dist.R")
# Pdf is saved: "coexpr_dist.pdf"

# Visualize co-expression network of HD-associated region in Gephi
source("gephiHDregion.R") #Requires large memory space

# Select genes co-expressed >0.5, for each polyQ gene and region
source("threshold_networks.R")
# Genes are also selectes when >0.4 or >0.6 to compare results.
# The number of genes based on each threshold is shown as a heatmap in "regionLs_threshold.pdf"

# Enrichment of GO terms for each set of genes that co-express with a polyQ gene in each region
source("goterms.R")

# Overlap of co-expressed genes between two polyQ genes


# Overlap of functional terms between genes co-expressed with polyQ genes



# Plot circular graphs of polyQ genes in Cytoscape, interactively
# source("cytPlot.R")






##### Function files #####
# tableNumbers.R
# volcanoplot.R
# setOverlap.R
# setOverlapSignif.R