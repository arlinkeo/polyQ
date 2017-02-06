# Number samples of atomic structures that fall within the HD associated regions (Coppen2016)

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/sampleIDs.RData")
#ontology <- read.csv("ABA_human_processed/Ontology_edited.csv")
load("resources/polyQ.RData")
# structureIDs[, 3] <- sapply(structureIDs[, 3], function(id){gsub(" ", "_", id)})
# structureIDs <- rbind(structureIDs, c(NA, "HDregion", "HD_region"))
# rownames(structureIDs) <- structureIDs$name
# probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
# entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element

# Number of samples in atomic structures and HD regions
donorNames <- names(sampleIDs[[1]])
sampleIDs_HD <- sampleIDs[["HD_region"]]
sampleIDs <- sampleIDs[-9]
overlap <- sapply(sampleIDs, function(s){
  sapply(donorNames, function(d){
    samples <- bitwAnd(s[[d]], sampleIDs_HD[[d]])
    #paste(length(which(samples == 1)), "/", length(which(s[[d]] == 1)), sep = "")
    #length(which(samples == 1)) / length(which(s[[d]] == 1)) * 100
    #length(which(samples == 1)) / length(which(sampleIDs_HD[[d]] == 1)) * 100
    length(which(samples == 1))
  })
})
overlap

# Number of samples in structure
sampleSize <- sapply(sampleIDs, function(s){
  length(which(unlist(s) == 1))
})
sampleSize