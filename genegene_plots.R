# Gene-gene expression scatterplots for all polyQ pairs across all donors.

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

load("resources/polyQ.RData")
load("resources/polyQ_expr.RData")

genepairs <- t(combn(polyQgenes, 2))
rownames(genepairs) <- apply(genepairs, 1, function(x){paste(x[1], "-", x[2], sep = "")})

lapply(names(donorList[1]), function(d){
  expr <- donorList[[d]]
  png(paste("genegene_plots_d", d, ".png", sep = ""), width = 1440, height = 1440, res = 150)
  par(mfrow=c(6,6))
  apply(genepairs, 1, function(x){
    gene1 <- expr[x[1], ]
    gene2 <- expr[x[2], ]
    plot(gene1, gene2)
  })
  dev.off()
})

