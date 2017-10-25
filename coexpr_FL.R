# co-expression matrix for frontal lobe

source("C:/Users/dkeo/surfdrive/polyQ_coexpression/PolyQ_scripts/baseScript.R")

load("resources/BrainExpr.RData")
load("resources/sampleIDs.RData")
corMatrix <- dget("PolyQ_scripts/corMatrix.R")
structureIDs <- rbind(HD_region = c(NA, "HDnetworkBD", "HD_region"), structureIDs)

#Specify brain region
id <- structureIDs["frontal_lobe", ]
sampleList <- sampleIDs[[id$name]]

#Select region-specific samples and calculate co-expression for each brain
brainCorList <- lapply(donorNames, function(d){
  samples <- sampleList[[d]] # row/col numbers to select
  expr <- brainExpr[[d]]
  corMatrix(samples, expr)
})
save(brainCorList, file = paste0("regional_coexpression/", structure, "/brainCorList_", id$acronym, ".RData"))
remove(brainExpr)
print("Correlation per brain saved")
names(brainCorList) <- NULL

#Correlation averaged across the six brains
meanCor <- apply(simplify2array(brainCorList), 1:2, mean)
save(meanCor, file = paste0("regional_coexpression/", structure, "/meanCor_", id$acronym, ".RData"))

print("Finished")
q(save = "no")