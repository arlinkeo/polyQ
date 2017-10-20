#Plotting region-specific networks
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression/regional_coexpression/")
library(WGCNA) # for heatmap function
options(stringsAsFactors = FALSE)

load("../resources/polyQ.RData")
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")
id <- ontology[ontology$name %in% "brain", ][ , c(1:3)] # change name for each region to be analyzed
#id <- c(0, "HDregion", "../HD_masks_Coppen2016")
load(paste(gsub(" ", "_", id[3]), "/subsetCor_", id[2], ".RData", sep = "")) # subset of top 25 correlated genes

#Heatmap of modules averaged across top 25 top correlated genes
modules <- rep(polyQgenes, each = 26)
moduleMat <- subsetCor
rownames(moduleMat) <- modules
colnames(moduleMat) <- modules
diag(moduleMat) <- NA

#Function to calculate mean correlation of pairwise modules, returns Module matrix with mean values
meanPerMod <- function (x) {
  mods <- unique(colnames(x))
  mat <- matrix(NA, length(mods), length(mods))
  rownames(mat) <- mods
  colnames(mat) <- mods
  for (r in mods) {
    for (c in mods) {
      modPair <- subsetCor[which(rownames(x) == r), which(colnames(x) == c)]
      mat[r, c] <- mean(modPair, na.rm = TRUE)
    }
  }
  mat
}

moduleMeans <- meanPerMod(moduleMat)
save(moduleMeans, file = paste(gsub(" ", "_", id[3]), "/moduleMeans_", id[2], ".RData", sep = ""))

make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}
labels <- as.vector(rbind(matrix("", 12, 9), polyQgenes, matrix("", 13, 9)))

for (i in 1:8) {
  pdf(file = paste(gsub(" ", "_", id[3]), "/polyQgenes_heatmap_", id[2], ".pdf", sep = ""), 8, 9)
  par(mar = c(5,6,12,3));
  labeledHeatmap(subsetCor, xLabels = make.italic(labels), colors = blueWhiteRed(200), zlim = c(-1,1), setStdMargins = FALSE, 
                 main = paste("polyQ top25 in ", id[3], sep =""), cex.lab = 1.3)

  labeledHeatmap(moduleMeans, xLabels = make.italic(colnames(moduleMeans)), colors = blueWhiteRed(200), zlim = c(-1,1), setStdMargins = FALSE, 
                 main = paste("Mean correlation within modules in ", id[3], sep =""), cex.lab = 1.3, textMatrix = round(moduleMeans, digits = 2))
  
  #Cluster modules
  geneTree = hclust(as.dist(1-moduleMeans), method = "average");
  order <- geneTree$order
  orderMods <- moduleMeans[order, order]

  par(mar = c(33.5,6.5,4,6.5), fig=c(0,1,0,1))
  plot(geneTree, xlab="", sub="", main = paste("Mean correlation of clustered modules in ", id[3], sep =""), 
       hang = 0.04, font = 3, axes = FALSE, ylab = "");
  par(mar = c(5,6,12,3), fig=c(0,1,0,1), new = TRUE)
  labeledHeatmap(orderMods, xLabels = make.italic(colnames(orderMods)), colors = blueWhiteRed(200), zlim = c(-1,1), 
                 setStdMargins = FALSE, textMatrix = round(orderMods, digits = 2))
  dev.off()
}