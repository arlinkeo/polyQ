# select correlated genes based on a corr. threshold

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellar nuclei","basal forebrain","globus pallidus"), ] # remove structures from list
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
  selection <- lapply(pQEntrezIDs, function(x){c(x, names(which(meanCor[x,] > 0.60)))})
  selection <- unlist(selection)
  matrix_selection <- meanCor[selection, selection]
  detach(2)
  matrix_selection
})
save(regionLs, file = "resources/regionLs_threshold060.RData")

regionLs <- lapply(structures, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/meanCor_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  selection <- lapply(pQEntrezIDs, function(x){c(x, names(which(meanCor[x,] > 0.70)))})
  selection <- unlist(selection)
  matrix_selection <- meanCor[selection, selection]
  detach(2)
  matrix_selection
})
save(regionLs, file = "resources/regionLs_threshold070.RData")

regionLs <- lapply(structures, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/meanCor_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  selection <- lapply(pQEntrezIDs, function(x){c(x, names(which(meanCor[x,] > 0.80)))})
  selection <- unlist(selection)
  matrix_selection <- meanCor[selection, selection]
  detach(2)
  matrix_selection
})
save(regionLs, file = "resources/regionLs_threshold080.RData")

##########################################################
load("regionLs_threshold080.RData")

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

make.table <- function(x){
  table <- sapply(x, function(y){
    pQidx <- which(colnames(y) %in% pQEntrezIDs)
    corr_genes <- sapply(c(1:9), function(i){# number of correlated genes
      if (i < 9){pQidx[i+1] - pQidx[i] - 1} #Minus the polyQ gene itself
      else {dim(y)[1] - pQidx[i]}
    })
    names(corr_genes) <- sapply(pQEntrezIDs, entrezId2Name)
    corr_genes
  })
  table <- rbind(table, apply(table, 2, sum))
  rownames(table)[10] <- "Total"
  table
}

roworder <- c("AR", "CACNA1A", "ATXN3", "ATN1", "HTT", "ATXN2", "ATXN1", "ATXN7", "TBP", "Total")
colorder <- c("striatum", "mesencephalon", "pons", "hypothalamus", "frontal_lobe", "cerebellar_cortex", "brain")

#Sort and plot table
pdf(file = "regionLs_threshold.pdf", 8, 9)
par(mar = c(2,6,12,3));
lapply(c(8, 7, 6), function(x){
  file <- paste("C:/Users/dkeo/surfdrive/polyQ_coexpression/resources/regionLs_threshold0", x, "0.RData", sep = "")
  attach(file)
  table <- make.table(regionLs)
  #sorted_table <- table[c(order(apply(table[-10, ], 1, sum), decreasing = TRUE), 10), order(table["Total", ], decreasing = TRUE)]
  sorted_table <- table[roworder, colorder]
  detach(2)
  labeledHeatmap(replace(sorted_table, which(sorted_table == 0), NA), xLabels = gsub("_", " ", colnames(sorted_table)), xLabelsPosition = "top", 
                 yLabels = c(make.italic(rownames(sorted_table)[-10]), rownames(sorted_table)[10]), colors = blueWhiteRed(200)[100:200], 
                 main = paste("Genes correlated >0.", x, " for each polyQ gene in different regions", sep = ), 
                 setStdMargins = FALSE, xLabelsAdj = 0, textMatrix = sorted_table)
})
dev.off()