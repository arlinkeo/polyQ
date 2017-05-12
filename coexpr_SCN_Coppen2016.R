# Co-expression analyses of brain regions (networks) from Coppen et al. 2016
# These regions correspond to structural co-variation networks based on grey matter differences from T1-MRI images between HD subjects and controls.

setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/polyQgenes/regional_coexpression/HD_region")
library(WGCNA)
allowWGCNAThreads(nThreads = 32)
options(stringsAsFactors = FALSE)

load("../../BrainExpr.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
load("../../sampleIDs.RData")
ontology <- read.csv("../../../../sjoerdhuisman/ABA_human_brain_probegene/Ontology_edited.csv")
load("../../polyQ.RData")
probeInfo <- read.csv("../../../../sjoerdhuisman/ABA_human_brain_probegene/probe_info_2014-11-11.csv")

#Select region-specific samples and calculate co-expression
brainCorList <- lapply(donorNames, function(d){
  samples <- as.logical(sampleIDs$HD_region[[d]]) # row/col numbers to select
  print(paste(d, ": ", sum(samples), " samples", sep = ""))
  expr <- brainExpr[[d]]
  corMat <- cor(t(expr[ , samples]))
  diag(corMat) <- 0
  corMat
})
save(brainCorList, file = "brainCorList_HDnetworkBD.RData")
remove(brainExpr)
print("Correlation per brain saved")
names(brainCorList) <- NULL

meanCor <- apply(simplify2array(brainCorList), 1:2, mean)
save(meanCor, file = "meanCor_HDnetworkBD.RData")
print("Finished")