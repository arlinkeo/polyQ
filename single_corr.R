# Single co-expression between polyQ genes in all regions incl. HD associated region (Coppen2016)
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
source("PolyQ_scripts/baseScript.R")
structureIDs <- rbind(c(NA, "HDnetworkBD", "HD_region"), structureIDs)
structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- structureIDs$name

### Load single correlations between polyQ genes ###
sc_list <- lapply(structures, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/meanCor_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  sc <- meanCor[pQEntrezIDs, pQEntrezIDs]
  colnames(sc) <- sapply(colnames(sc), entrezId2Name)
  rownames(sc) <- sapply(rownames(sc), entrezId2Name)
  detach(2)
  sc
})
save(sc_list, file = "resources/sc_list.RData")

###PLot heatmaps for all regions
structureIDs <- structureIDs[!structureIDs$name %in% c("brain", "cerebellum"), ]
mat_hd <- sc_list$HD_region#Cluster in HD
geneTree = hclust(as.dist(1-mat_hd), method = "average");
order <- rev(geneTree$order)#Order of all maps as in HD region

pdf(file = "coexpression_heatmap.pdf", 8, 9)
par(mar = c(5,6,12,3));
lapply(structureIDs$name, function(r){
  mat <- sc_list[[r]][order, order]
  labels <- rownames(mat)
  labeledHeatmap(mat, xLabels = labels, colors = blueWhiteRed(200), zlim = c(-1,1), setStdMargins = FALSE, 
               main = paste("Direct co-expression in ", r, sep =""), cex.lab = 1.3, textMatrix = round(mat, digits = 2))
})
dev.off()