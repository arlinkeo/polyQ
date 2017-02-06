##### Obtain binary vectors indicating presence of sample across regions and donors #####
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/BrainExpr.RData")
ontology <- read.csv("ABA_human_processed/Ontology_edited.csv")
load("resources/polyQ.RData")
structureIDs[, 3] <- sapply(structureIDs[, 3], function(id){gsub(" ", "_", id)})
structureIDs <- rbind(structureIDs, c(NA, "HDregion", "HD_region"))
rownames(structureIDs) <- structureIDs$name
donorNames <- names(brainExpr)
names(donorNames) <- donorNames

#Select HD region-specific samples
sampleIDs_HD <- lapply(donorNames, function(d){
  networkInfo <- read.csv(paste("regional_coexpression/HD_region/samples_in_networks_", d, ".txt", sep = ""), header = TRUE, sep = "\t", comment.char = "#")
  networkB <- networkInfo[, "network_B"]
  networkD <- networkInfo[, "network_D"]
  bitwOr(networkB, networkD)
})
#Select anatomic region-specific samples
sampleIDs <- apply(structureIDs[-9, ], 1, function(id){
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
#Combine info
sampleIDs <- c(sampleIDs, HD_region = list(sampleIDs_HD))
save(sampleIDs, file = "resources/sampleIDs.RData")