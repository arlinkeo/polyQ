# Number samples of atomic structures that fall within the HD associated regions (Coppen2016)
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)
library(ggplot2)

# load("resources/polyQ.RData")
# structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellum"), ]
#structureIDs <- rbind(HD_region = c(NA, "HDregion", "HD_region"), structureIDs)
load("resources/brainExpr.RData") 
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
load("resources/sampleIDs.RData")
sampleIDs_HD <- sampleIDs[["HD_region"]]
sampleIDs <- sampleIDs[!names(sampleIDs) %in% c("HD_region", "brain", "cerebellum")]


# Number of samples in atomic structures and HD regions
overlap <- lapply(sampleIDs, function(s){
  lapply(donorNames, function(d){
    samples <- bitwAnd(s[[d]], sampleIDs_HD[[d]])
  })
})
tab <- sapply(overlap, function(d){sapply(d, sum)})
totalN <- apply(tab, 2, sum)
# 
# Total samples in and out HD region
totalHD <- sum(sapply(sampleIDs_HD, sum))
idHD <- unlist(sampleIDs_HD)
notHD <- as.numeric(!idHD)

allStruct <- lapply(sampleIDs, unlist)
allStruct <-Reduce(bitwOr, allStruct)
nonStruct <- as.numeric(!allStruct)
restRest <- bitwAnd(nonStruct, notHD)

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

recur <- function(x){
  if (x == 1 ) tab$samples[x]
  else tab$samples[x]+recur(x-1)
}
tab$ymax <- sapply(c(1:nrow(tab)), recur)
tab$ymin <- tab$ymax-tab$samples


ggplot(tab) + geom_rect(aes(fill = HD), ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)
+ geom_rect(aes(fill = structure), ymax = ymax, ymin = ymin, xmax = 3, xmin = )
# # Number of samples in Hd region minus structure, should add up to total number samples in structures
# diffHD <- lapply(sampleIDs[-1], function(s){
#   lapply(donorNames, function(d){
#     overlapping <- bitwAnd(s[[d]], sampleIDs_HD[[d]])
#     samples <- sampleIDs_HD[[d]] - overlapping
#   })
# })
# sapply(diffHD, function(d){sapply(d, sum)})
# # 
# apply(structureIDs, 1, function(r){
#   structure <- r$acronym
#   structName <- r$name
#   overlapSamples <- overlap[[structName]]
#   diffHDsamples <- diffHD[[structName]]
#   
#   brainCorList <- lapply(donorNames, function(d){
#     expr <- brainExpr[[d]][, overlapSamples[[d]]]
#     corMat <- cor(t(expr))
#     diag(corMat) <- 0
#     corMat
#   })
#   
# })