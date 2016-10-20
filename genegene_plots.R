# Gene-gene expression scatterplots for all polyQ pairs across all donors.

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

load("resources/polyQ.RData")
ontology <- read.csv("ABA_human_processed/Ontology_edited.csv")
rownames(ontology) <- ontology$id
load("resources/polyQ_expr.RData")

genepairs <- t(combn(polyQgenes, 2))
rownames(genepairs) <- apply(genepairs, 1, function(x){paste(x[1], "-", x[2], sep = "")})

lapply(names(donorList[1]), function(d){
  expr <- donorList[[d]]
  colors <- sapply(colnames(expr), function(x){
    clr <- ontology[x, 'color_hex_triplet']
    if (nchar(clr) == 5) {paste("#0", clr, sep = "")}
    else {paste("#", clr, sep = "")}
  })
  png(file = paste("genegene_plots_d", d, ".png", sep = ""), width = 1920, height = 1920, res = 300, pointsize = 4)
  par(mfrow=c(6,6))
  apply(genepairs, 1, function(x){
    gene1 <- expr[x[1], ]
    gene2 <- expr[x[2], ]
    plot(gene1, gene2, xlab = x[1], ylab = x[2], col = colors, cex = 0.5, cex.lab = 2)
  })
  dev.off()
})

