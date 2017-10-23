# Number of samples of atomic structures that fall within the HD associated region (Coppen2016)

source("C:/Users/dkeo/surfdrive/polyQ_coexpression/PolyQ_scripts/baseScript.R")

load("resources/brainExpr.RData") 
load("resources/sampleIDs.RData")
sampleIDs_HD <- sampleIDs[["HD_region"]]
sampleIDs <- sampleIDs[!names(sampleIDs) %in% c("HD_region", "brain")]

# Number of samples in atomic structures and HD regions
overlap <- lapply(sampleIDs, function(s){
  lapply(donorNames, function(d){
    samples <- bitwAnd(s[[d]], sampleIDs_HD[[d]])
  })
})
tab <- sapply(overlap, function(d){sapply(d, sum)})
tab
totalN <- apply(tab, 2, sum)
totalN

# Total samples in and out HD region
totalHD <- sum(sapply(sampleIDs_HD, sum)) # total number of samples in HD-associated region
idHD <- unlist(sampleIDs_HD)
notHD <- as.numeric(!idHD)

# Total samples in anatomical regions
allStruct <- lapply(sampleIDs, unlist)
allStruct <-Reduce(bitwOr, allStruct)
nonStruct <- as.numeric(!allStruct)
restRest <- bitwAnd(nonStruct, notHD) # Samples not in anatomical structure and not in HD-associated region

# For each anatomical region, number of samples in- and outside HD-associated region
colHD <- lapply(names(sampleIDs), function(n){
  s <- sampleIDs[[n]]
  inHD <- sum(bitwAnd(unlist(s), idHD))
  nonHD <- sum(bitwAnd(unlist(s), as.numeric(!idHD)))
  inHD <- c("HD-associated region", n, inHD)
  nonHD <- c("rest", n, nonHD)
  rbind(inHD, nonHD)
})
tab <- Reduce(rbind, colHD)
rownames(tab) <- NULL
tab <- gsub("_", " ", tab)
tab <- as.data.frame(tab)
colnames(tab) <- c("HD", "structure", "samples")
tab[nrow(tab)+1, ] <- c("rest", "rest", sum(restRest))
tab$samples <- as.numeric(tab$samples)
tab

# Number of samples in Hd region minus structure, should add up to total number samples in HD-region
diffHD <- lapply(sampleIDs, function(s){
  lapply(donorNames, function(d){
    overlapping <- bitwAnd(s[[d]], sampleIDs_HD[[d]])
    samples <- sampleIDs_HD[[d]] - overlapping
  })
})
sapply(diffHD, function(d){sapply(d, sum)})