# select correlated genes based on a corr. threshold

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

#####Load Data#########

#Load mean corr. data across 6 brains. and select based on threshold
structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- gsub(" ", "_", structureIDs$name)

regionLs <- lapply(structures, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/meanCor_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  selection <- lapply(pQEntrezIDs, function(x){c(x, names(which(meanCor[x,] > 0.50)))})
  names(selection) <- pQEntrezIDs
  detach(2)
  selection
})
save(regionLs, file = "resources/genesets_threshold050.RData")

regionLs <- lapply(structures, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/meanCor_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  selection <- lapply(pQEntrezIDs, function(x){c(x, names(which(meanCor[x,] > 0.60)))})
  names(selection) <- pQEntrezIDs
  detach(2)
  selection
})
save(regionLs, file = "resources/genesets_threshold060.RData")

regionLs <- lapply(structures, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/meanCor_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  selection <- lapply(pQEntrezIDs, function(x){c(x, names(which(meanCor[x,] > 0.70)))})
  names(selection) <- pQEntrezIDs
  detach(2)
  selection
})
save(regionLs, file = "resources/genesets_threshold070.RData")

regionLs <- lapply(structures, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/meanCor_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  selection <- lapply(pQEntrezIDs, function(x){c(x, names(which(meanCor[x,] > 0.80)))})
  names(selection) <- pQEntrezIDs
  detach(2)
  selection
})
save(regionLs, file = "resources/genesets_threshold080.RData")

##########################################################
load("resources/genesets_threshold050.RData")

#Export network to Cytoscape and plot heatmaps
lapply(c(1:9), function(x){
  name <- allID[x, ]
  mat <- regionLs[[x]]
  rownames(mat) <- sapply(rownames(mat), entrezId2Name)
  colnames(mat) <- sapply(colnames(mat), entrezId2Name)
  pQidx <- which(colnames(mat) %in% polyQgenes)
  colors <- sapply(1:dim(mat)[1], function(i){
    if (i %in% pQidx) {pQcolors[which(pQidx == i)]}
    else { while(!(i %in% pQidx)) {i <- i-1}; pQcolors[which(pQidx == i)]}
  })
  names(colors) <- colnames(mat)
  cyt <- exportNetworkToCytoscape(mat, threshold = 0.0, edgeFile = paste(name[3], "/Threshold080Edges_", name[2],".txt", sep = ""),
                                  nodeFile = paste(name[3], "/Threshold080Nodes_", name[2],".txt", sep = ""), nodeAttr = colors)
  
  # pdf(file = paste(name[3], "/Threshold080_heatmap_", gsub("_", " ", name[2]), ".pdf", sep = ""), 8, 9)
  # par(mar = c(5,6,12,3));
  # rownames(mat) <- sapply(rownames(mat), entrezId2Name)
  # rn_pq <- which(colnames(mat) %in% polyQgenes)
  # labels <- replace(colnames(mat), -rn_pq, "")
  # labeledHeatmap(mat, xLabels = labels, colors = blueWhiteRed(200), setStdMargins = FALSE, 
  #                main = paste("Genes correlated >0.8 with a polyQ gene in the", name[3]), zlim = c(-1,1))
  # dev.off()
})

#Overlapping gene-sets
dups <- sapply(regionLs, function(x){idx <- which(duplicated(colnames(x))); sapply(colnames(x)[idx], entrezId2Name)}) 

#sum <- sapply(regionLs, function(x){dim(x)[1]-9}) #Minus 9 polyQ genes

######################################################################################################
#Number of genes correlated for each polyQ gene in different regions with different thresholds
#roworder <- c("ATN1", "CACNA1A", "HTT", "AR", "ATXN2", "ATXN1", "ATXN3", "ATXN7", "TBP", "Total")
colorder <- c("brain", "frontal_lobe", "parietal_lobe", "striatum", "hypothalamus", "mesencephalon", "cerebellar_cortex", "pons")

#Sort and plot table
pdf(file = "regionLs_threshold.pdf", 8, 9)
par(mar = c(2,6,12,3));
lapply(c(5:8), function(x){
  file <- paste("C:/Users/dkeo/surfdrive/polyQ_coexpression/resources/genesets_threshold0", x, "0.RData", sep = "")
  attach(file)
  table <- sapply(regionLs, function(x){sapply(x, length)})
  detach(2)
  rownames(table) <- sapply(rownames(table), entrezId2Name)
  Total <- apply(table, 2, sum)
  table <- rbind(table, Total)
  table <- table[ , colorder]
  labeledHeatmap(replace(table, which(table == 1), NA), xLabels = gsub("_", " ", colnames(table)), xLabelsPosition = "top", 
                 yLabels = c(make.italic(rownames(table)[-10]), rownames(table)[10]), colors = blueWhiteRed(200)[100:200], 
                 main = paste("Genes correlated >0.", x, " for each polyQ gene in different regions", sep = ""), 
                 setStdMargins = FALSE, xLabelsAdj = 0, textMatrix = table)
})
dev.off()