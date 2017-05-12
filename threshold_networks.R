# select correlated genes based on a corr. threshold

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("brain", "cerebellum"), ]
structureIDs <- rbind(HD_region = c(NA, "HDnetworkBD", "HD_region"), structureIDs)
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
# make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

#Load mean corr. data across 6 brains. and select based on threshold
structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- structureIDs$name

# for(i in seq(0.4, 0.6, by = 0.1)){
#   regionLs <- lapply(structures, function(x) {
#     f <- paste("regional_coexpression/", x[3], "/meanCor_", x[2], ".RData", sep = "")
#     print(f)
#     attach(f)
#     selection <- lapply(pQEntrezIDs, function(x){c(x, names(which(meanCor[x,] > i)))})
#     names(selection) <- pQEntrezIDs
#     detach(2)
#     selection
#   })
#   save(regionLs, file = paste("resources/genesets_threshold0", i*10,"0.RData", sep = ""))
# }

#Export network to Cytoscape and plot heatmaps
# lapply(c(1:9), function(x){
#   name <- allID[x, ]
#   mat <- regionLs[[x]]
#   rownames(mat) <- sapply(rownames(mat), entrezId2Name)
#   colnames(mat) <- sapply(colnames(mat), entrezId2Name)
#   pQidx <- which(colnames(mat) %in% polyQgenes)
#   colors <- sapply(1:dim(mat)[1], function(i){
#     if (i %in% pQidx) {pQcolors[which(pQidx == i)]}
#     else { while(!(i %in% pQidx)) {i <- i-1}; pQcolors[which(pQidx == i)]}
#   })
#   names(colors) <- colnames(mat)
#   cyt <- exportNetworkToCytoscape(mat, threshold = 0.0, edgeFile = paste(name[3], "/Threshold080Edges_", name[2],".txt", sep = ""),
#                                   nodeFile = paste(name[3], "/Threshold080Nodes_", name[2],".txt", sep = ""), nodeAttr = colors)
#   
#   # pdf(file = paste(name[3], "/Threshold080_heatmap_", gsub("_", " ", name[2]), ".pdf", sep = ""), 8, 9)
#   # par(mar = c(5,6,12,3));
#   # rownames(mat) <- sapply(rownames(mat), entrezId2Name)
#   # rn_pq <- which(colnames(mat) %in% polyQgenes)
#   # labels <- replace(colnames(mat), -rn_pq, "")
#   # labeledHeatmap(mat, xLabels = labels, colors = blueWhiteRed(200), setStdMargins = FALSE, 
#   #                main = paste("Genes correlated >0.8 with a polyQ gene in the", name[3]), zlim = c(-1,1))
#   # dev.off()
# })

#Sort and plot table
colorder <- c("HD_region", "frontal_lobe", "parietal_lobe", "striatum", "hypothalamus", "mesencephalon", "cerebellar_cortex", "pons")

pdf(file = "regionLs_threshold.pdf", 12, 4.5)
par(mfrow = c(1, 3), oma = c(1, 2, 0, 1), lwd = 0.01, mai = c(0, .6, 1.5, .3));
for(i in 4:6){
  f <- paste("resources/genesets_threshold0", i, "0.RData", sep = "")
  attach(f)
  table <- sapply(regionLs, function(x){sapply(x, length)})
  rownames(table) <- sapply(rownames(table), entrezId2Name)
  detach(2)
  labeledHeatmap(replace(table, which(table == 1), NA), xLabels = gsub("_", " ", colnames(table)), xLabelsPosition = "top", 
               yLabels = rownames(table), colors = blueWhiteRed(200)[100:200], 
               main = paste("co-expression >0.", i, "0", sep = ""), 
               setStdMargins = FALSE, xLabelsAdj = 0, textMatrix = table, cex.main = 1.2, cex.text = 1, cex.lab = 1)
  }
dev.off()