# Number samples of atomic structures that fall within the HD associated regions (Coppen2016)

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/sampleIDs.RData")
load("resources/polyQ.RData")

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