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

#Select region-specific samples
sample_ids <- lapply(donorNames, function(d){
  networkInfo <- read.csv(paste("samples_in_networks_", d, ".txt", sep = ""), header = TRUE, sep = "\t", comment.char = "#")
  networkBSamples <- networkInfo[which(networkInfo[, "network_B"] == 1), "sample_id"]
  networkDSamples <- networkInfo[which(networkInfo[, "network_D"] == 1), "sample_id"]
  combinedSamples <- as.character(unique(c(networkBSamples, networkDSamples)))
})
lapply(sample_ids, length)

brainCorList <- lapply(donorNames, function(d) {
  samples <- sample_ids[[d]]
  expr <- brainExpr[[d]]
  cols <- colnames(expr) %in% samplesprint(paste("Selected samples: ", length(which(cols))))
  corMat <- cor(t(expr[ , cols]))
  diag(corMat) <- 0
  corMat
})
save(brainCorList, file = paste("brainCorList_HDnetworkBD.RData", sep = ""))
remove(brainExpr)
print("Correlation per brain saved")
load("../brainCorList_HDnetworkBD.RData")
names(brainCorList) <- NULL

meanCor <- apply(simplify2array(brainCorList), 1:2, mean)
save(meanCor, file = "meanCor_HDnetworkBD.RData")
print("Mean corr. accross brains saved")

print("Finished")