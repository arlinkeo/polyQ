# Number samples of atomic structures that fall within the HD associated regions (Coppen2016)

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
#library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/BrainExpr.RData")
ontology <- read.csv("ABA_human_processed/Ontology_edited.csv")
load("resources/polyQ.RData")
rownames(structureIDs) <- structureIDs$name
donorNames <- names(brainExpr)
names(donorNames) <- donorNames

#Select HD region-specific samples
sampleIDs_HD <- lapply(donorNames, function(d){
  networkInfo <- read.csv(paste("HD_masks_Coppen2016/samples_in_networks_", d, ".txt", sep = ""), header = TRUE, sep = "\t", comment.char = "#")
  networkB <- networkInfo[, "network_B"]
  networkD <- networkInfo[, "network_D"]
  bitwOr(networkB, networkD)
  # networkBSamples <- which(networkInfo[, "network_B"] == 1)
  # networkDSamples <- which(networkInfo[, "network_D"] == 1)
  # combinedSamples <- sort(unique(c(networkBSamples, networkDSamples))) # row/col numbers to select
  # print(paste(d, ": ", length(combinedSamples), " samples", sep = ""))
  # combinedSamples
})

#Select anatomic region-specific samples
sampleIDs <- apply(structureIDs, 1, function(id){
  print(id)
  structName <- id[2]
  ontologyRows <- grep(id[1], ontology$structure_id_path)
  selectIds <- as.character(ontology$id[ontologyRows])
  lapply(donorNames, function(d){
    expr <- brainExpr[[d]]
    ids <- intersect(selectIds, colnames(expr))
    cols <- colnames(expr) %in% ids
    print(paste("Samples: ", length(which(cols))))
    as.numeric(cols)
  })
})

# Number of samples in atomic structures and HD regions
overlap <- sapply(sampleIDs, function(s){
  sapply(donorNames, function(d){
    samples <- bitwAnd(s[[d]], sampleIDs_HD[[d]])
    #paste(length(which(samples == 1)), "/", length(which(s[[d]] == 1)), sep = "")
    #length(which(samples == 1)) / length(which(s[[d]] == 1)) * 100
    #length(which(samples == 1)) / length(which(sampleIDs_HD[[d]] == 1)) * 100
    length(which(samples == 1))
  })
})
