# Number samples of atomic structures that fall within the HD associated regions (Coppen2016)
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

# load("resources/polyQ.RData")
# structureIDs <- structureIDs[!structureIDs$name %in% c("brain", "cerebellum"), ]
#structureIDs <- rbind(HD_region = c(NA, "HDregion", "HD_region"), structureIDs)
load("resources/brainExpr.RData") 
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
load("resources/sampleIDs.RData")
sampleIDs <- sampleIDs[!names(sampleIDs) %in% c("cerebellum")]
sampleIDs_HD <- sampleIDs[["HD_region"]]

# Number of samples in atomic structures and HD regions
overlap <- lapply(sampleIDs[-1], function(s){
  lapply(donorNames, function(d){
    samples <- bitwAnd(s[[d]], sampleIDs_HD[[d]])
  })
})
sapply(overlap, function(d){sapply(d, sum)})
# 
# # Number of samples in Hd region minus structure, should add up to toal number samples in structures
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