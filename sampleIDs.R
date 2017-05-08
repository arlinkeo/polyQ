##### Obtain binary vectors indicating presence of sample across regions and donors #####
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/BrainExpr.RData")
ontology <- read.csv("ABA_human_processed/Ontology_edited.csv")
load("resources/polyQ.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames

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
    as.integer(cols)
  })
})
#Select HD region-specific samples
structureIDs <- rbind(structureIDs, c(NA, "HDregion", "HD_region"))
sampleIDs_HD <- lapply(donorNames, function(d){
  networkInfo <- read.csv(paste("HD_mask/Huntington_results/networksamples_", d, ".txt", sep = ""), header = TRUE, sep = "\t", comment.char = "#")
  as.numeric(networkInfo$Inside.y.n)
})
#Combine info
sampleIDs <- c(HD_region = list(sampleIDs_HD), sampleIDs)
save(sampleIDs, file = "resources/sampleIDs.RData")