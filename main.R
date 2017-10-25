# Main script
# All scripts are in directory PolyQ_scripts

# Obtain data information of polyQ genes and brain regions
source("PolyQ_scripts/polyQ_region.R")
# Information is stored in "resources/polyQ.RData" and called whenever needed

# baseScript.R is called in every script to have basic function an variables needed

# Read microarray expression data and save as R-object 
source("PolyQ_scripts/brainExpr.R")
# Data saved as "resources/BrainExpr.RData"

# Select region-specific samples
# Save as binary vectors indicating for each sample whether it is within the specified region (1) or not (0)
source("PolyQ_scripts/sampleIDs.R")
# Sample selection saved as "resources/sampleIDs.RData"

# Plot pie charts of samples in HD region and samples in anatomical regions
source("PolyQ_scripts/pieChart_samples.R")

# Check number of samples of atomic structures that fall within the HD associated region (Coppen2016)
source("PolyQ_scripts/overlap_samples.R")

# Mean and variance of polyQ gene expression across regional-specific samples and all donors
source("PolyQ_scripts/avgExpr.R")
# Mean expr. is saved as "resources/avgExpr.RData"
# Mean and variance table is saved as heatmap in "expr_meanvarmedian.pdf"
# Color is assigned to the mean expression and saved as "resources/avgExprColor.RData"
# Legend bar of the expression is saved in "avgExpr.pdf"

# Differential expression of HD-associated region
# Analysis is done for each brain separately, then genes are selected if significant in at least 5 brains
# For the gene selection, functional enrichment analysis is doen for up- and down-regulated genes
source("PolyQ_scripts/diffExpr_HDregion.R")
# Differentially expressed genes are saved as "resources/diffgenesPerBrain.RData"
# A pdf is saved with  boxplots of the expression of each polyQ genes in- and outside the HD-region.

# Co-expression network of HD-associated region and anatomical brain regions
# Separate script for each brain region, because of memory requirements
source("PolyQ_scripts/coexpr_HDregion.R")
source("PolyQ_scripts/coexpr_FL.R")
source("PolyQ_scripts/coexpr_PL.R")
source("PolyQ_scripts/coexpr_Str.R")
source("PolyQ_scripts/coexpr_Hy.R")
source("PolyQ_scripts/coexpr_MES.R")
source("PolyQ_scripts/coexpr_CbCx.R")
# correlation matrices saved in separate directories for each brain region
# "brainCorList_[region].RData" and "meanCor_[region].RData"

# Extract direct, single correlation between polyQ genes in each region
source("PolyQ_scripts/single_corr.R")
# Saved as "resources/sc_list.RData"

# Plot histograms of correlation distribution of each polyQ gene and brain region
source("PolyQ_scripts/coexpr_dist.R")
# Pdf is saved: "coexpr_dist.pdf"

# Visualize co-expression network of HD-associated region in Gephi
source("PolyQ_scripts/gephiHDregion.R") #Requires large memory space

# Select genes co-expressed >0.5, for each polyQ gene and region
source("PolyQ_scripts/threshold_networks.R")
# Genes are also selectes when >0.4 or >0.6 to compare results.
# The number of genes based on each threshold is shown as a heatmap in "regionLs_threshold.pdf"

# Enrichment of GO terms for each set of genes that co-express with a polyQ gene in each region
source("PolyQ_scripts/goterms.R")
# results saved in folder for each brain region: regional_coexpression/[brain region]/goterms050_[brain region acronym]_[polyQ gene].txt

# Overlap of co-expressed genes between two polyQ genes
source("PolyQ_scripts/overlap_genesets.R")
# saved as "resources/geneSetOverlap.RData" and "resources/geneSetOverlapSignif.RData"
# Furthermore, the presence of dna repair and ubiquitination genes is  checked

# Overlap of functional terms between genes co-expressed with polyQ genes
source("PolyQ_scripts/overlap_goterms.R")
# Number of terms found for each brain region and polyQ gene is saved as "number_of_goterms050.pdf"
# Terms that overlap between polyQ gene sets are saved as "resources/termSetOverlap.RData"
# those that overlap with HTT, ATN1, ATXN2 in the HD-associated region are saved as a table: "overlapTerms_HTT_ATN1_ATXN2.txt"

# Plot overlap of genes and terms for each gene pair
source("PolyQ_scripts/plot_overlap_table.R")

# Plot circular graphs of polyQ genes in Cytoscape, interactively
# source("PolyQ_scripts/cytPlot.R")






##### Function files #####
# tableNumbers.R
# volcanoplot.R
# setOverlap.R
# setOverlapSignif.R