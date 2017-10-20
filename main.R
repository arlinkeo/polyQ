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



