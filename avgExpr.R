#Average expression of polyQ gene across regional-specific samples and all donors

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)
library(RColorBrewer)
load("resources/polyQ.RData")
# load("resources/BrainExpr.RData")
# load("resources/sampleIDs.RData")
# sampleIDs$brain <- NULL
# sampleIDs$cerebellum <- NULL
# donorNames <- names(brainExpr)
# names(donorNames) <- donorNames
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element

#Average expression of a polyQ gene in a structure across donors and samples.
# avgExpr <- sapply(sampleIDs, function(s){
#   res <- sapply(donorNames, function(d){
#     expr <- brainExpr[[d]]
#     expr2 <- expr[pQEntrezIDs, as.logical(s[[d]])]
#     apply(expr2, 1, mean) # Avg across region-specific samples per PQ per donor
#   })
#   apply(res, 1, mean) # Avg across region-specific samples and donors per PQ
# })
# rownames(avgExpr) <- sapply(rownames(avgExpr), entrezId2Name)
# save(avgExpr, file = "resources/avgExpr.RData")
load("resources/avgExpr.RData")

pal <- colorRampPalette(rev(brewer.pal(3, "RdBu")))# pal = colorRampPalette(c('darkgreen', 'yellow', 'red'))
avgExprColor <- matrix(pal(100)[as.numeric(cut(avgExpr, breaks = 100))], length(polyQgenes), ncol(avgExpr))
rownames(avgExprColor) <- rownames(avgExpr)
colnames(avgExprColor) <- colnames(avgExpr)
save(avgExprColor, file = "resources/avgExprColor.RData")