# Co-expression analyses of brain regions (networks) from Coppen et al. 2016
# These regions correspond to structural co-variation networks based on grey matter differences from T1-MRI images between HD subjects and controls.

setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/polyQgenes/HD_masks_Coppen2016")
library(WGCNA)
allowWGCNAThreads(nThreads = 32)
options(stringsAsFactors = FALSE)

load("../BrainExpr.RData")
donorNames <- list("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

ontology <- read.csv("../../../sjoerdhuisman/ABA_human_brain_probegene/Ontology_edited.csv")
load("../polyQ.RData")
probeInfo <- read.csv("../../../sjoerdhuisman/ABA_human_brain_probegene/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]}

#Select region-specific samples and calculate co-expression
brainCorList <- lapply(donorNames, function(d){
  networkInfo <- read.csv(paste("samples_in_networks_", d, ".txt", sep = ""), header = TRUE, sep = "\t", comment.char = "#")
  networkBSamples <- which(networkInfo[, "network_B"] == 1)
  networkDSamples <- which(networkInfo[, "network_D"] == 1)
  combinedSamples <- sort(unique(c(networkBSamples, networkDSamples))) # row/col numbers to select
  print(paste(d, ": ", length(combinedSamples), " samples", sep = ""))
  expr <- brainExpr[[d]]
  corMat <- cor(t(expr[ , combinedSamples]))
  diag(corMat) <- 0
  corMat
})
save(brainCorList, file = paste("brainCorList_HDnetworkBD2.RData", sep = ""))
remove(brainExpr)
print("Correlation per brain saved")
#load("brainCorList_HDnetworkBD2.RData")
names(brainCorList) <- NULL

meanCor <- apply(simplify2array(brainCorList), 1:2, mean)
save(meanCor, file = "meanCor_HDnetworkBD2 .RData")
print("Mean corr. accross brains saved")

load("../HD_masks_Coppen2016/meanCor_HDnetworkBD.RData")
top25id <- apply(meanCor[pQEntrezIDs, ], 1, function(x) {names(head(-sort(-x), 25))})
top25id <- rbind(colnames(top25id), top25id)
top25sym <- apply(top25id,c(1,2), entrezId2Name)
geneSet <- as.vector(top25id) # May contain duplicates
subsetCor <- meanCor[geneSet, geneSet]
geneSetSym <- sapply(geneSet, entrezId2Name)
rownames(subsetCor) <- geneSetSym
colnames(subsetCor) <- geneSetSym
save(subsetCor, file = "subsetCor_HDregion.RData")

print("Finished")