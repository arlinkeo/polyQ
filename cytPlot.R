# Plot data in Cytoscape
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(RCy3)
library(reshape2)
options(stringsAsFactors = FALSE)

source("PolyQ_scripts/baseScript.R")
load("resources/avgExpr.RData")
load("resources/avgExprColor.RData")

structureIDs <- rbind(c(NA, "HDregion", "HD_region"), structureIDs)
rownames(structureIDs) <- structureIDs$name

pqOrder <- c(3, c(1:2), c(4:9))
polyQgenes <- polyQgenes[pqOrder]
pQEntrezIDs <- pQEntrezIDs[pqOrder]
pQcolors <- pQcolors[pqOrder]
avgExpr <- avgExpr[pqOrder, ]
avgExprColor <- avgExprColor[pqOrder, ]

pqPairs <- t(combn(polyQgenes, 2))
rownames(pqPairs) <- apply(pqPairs, 1, function(x){paste(x[1], "-", x[2], sep = "")})

#Functions
linMap <- function(x){
  if (x>0.5) round((x - 0.5) * 2 * 10 + 2, digits = 1)
  else 0
}
# Cytoscape default visuals
cyt.visuals <- function(cw) {
  layoutNetwork (cw, layout.name = "circular")
  setDefaultNodeShape (cw, "Ellipse")
  setDefaultNodeFontSize (cw, 16)
  #lockNodeDimensions (cw, TRUE) # Not implemented yet in RCy3
  #setDefaultNodeSize(cw, 35)
  #showGraphicsDetails(cw, TRUE) # Not implemented yet in RCy3
}

#re-order function for table with polyQ gene pairs as rows
new.order <- function(o, data){
  pairs <- t(sapply(rownames(data), function(pair){unlist(strsplit(pair, split = "-"))}))
  res <- t(apply(o, 1, function(pg){
    row <- which(apply(pairs, 1, function(g){
      Reduce("&", pg %in% g)
    }))
    t(data[row, ])
  }))
  colnames(res) <- colnames(data)
  res
}

# Load info about direct interacters (single_corr.R)
load("resources/sc_list.RData")
#min expr across regions
min(sapply(sc_list[names(sc_list) != "cerebellum"], function(x){x[x<=0.5]= NA;min(x, na.rm = TRUE)}))
# Plot pqNeighbors in Cytoscape for each structure
apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  mat <- sc_list[[structure]][pqOrder, pqOrder] # adjacency matrix
  mat[lower.tri(mat)] <- 0
  edgeTable <- as.data.frame(as.table(mat))
  colnames(edgeTable) <- c("fromNode", "toNode", "coexpr")
  edgeTable$edgeType <- rep("coexpr", dim(edgeTable)[1])
  edgeTable$width <- sapply(edgeTable$coexpr, linMap)
  nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], color = avgExprColor[, structure])
  # nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], color = pQcolors)

  g <- cyPlot(nodeTable, edgeTable)
  windowName <- paste("pqNeighbors_",structName, sep = "")
  cw <- CytoscapeWindow(windowName, graph = g, overwrite = TRUE)
  displayGraph (cw)
  cyt.visuals(cw)
  setNodeColorDirect(cw,  nodeTable[ , "nodeName"],  nodeTable[ , "color"])
  edgeRows <- which(edgeTable$coexpr > 0.5)
  edges <- apply(edgeTable[edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (coexpr) ")})
  nonEdges <- apply(edgeTable[-edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (coexpr) ")})
  selectEdges(cw, nonEdges)
  deleteSelectedEdges(cw)
  setEdgeLineWidthDirect(cw, edges, edgeTable$width[edgeRows])
  # svgName <- paste("Images/pqNeighbors2_", structure, "/", windowName, ".svg", sep = "")
  # saveImage(cw, svgName, "svg")
})

# Load new cytoscape session

# Gene set overlap info
load("resources/geneSetOverlap.RData")
overlapNumbers <- data.frame(sapply(geneSetOverlap, function(r){sapply(r, length)}))
overlapNumbers <- new.order(pqPairs, overlapNumbers)
load("resources/geneSetOverlapSignif.RData")
geneSetOverlapSignif <- new.order(pqPairs, geneSetOverlapSignif)
#max and min across regions
overlapNumbers[geneSetOverlapSignif >= 0.05] <- NA
max(apply(overlapNumbers, 2, max))
min(apply(overlapNumbers, 2, function(x) min(x, na.rm = TRUE)));
# Plot for overlap of gene sets
interaction_list <- apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  print(structure)
  overlap <- as.numeric(overlapNumbers[, structure])
  width <- sapply(overlap, log1p)
  signif <- geneSetOverlapSignif[, structure] # adjacency matrix
  
  nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], color = avgExprColor[, structure])
  # nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], color = pQcolors)
  edgeTable <- data.frame(fromNode = pqPairs[ , 1], toNode = pqPairs[ , 2], edgeType="overlap",
                          overlap = overlap, width = width, signif = signif)
  g <- cyPlot(nodeTable, edgeTable)
  windowName <- paste("overlapGeneSets_",structName, sep = "")
  cw <- CytoscapeWindow(windowName, graph = g, overwrite = TRUE)
  displayGraph (cw)
  cyt.visuals(cw)
  setNodeColorDirect(cw,  nodeTable[ , "nodeName"],  nodeTable[ , "color"])
  
  edgeRows <- which(edgeTable$signif < 0.05)
  edges <- apply(edgeTable[edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (overlap) ")})
  nonEdges <- apply(edgeTable[-edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (overlap) ")})
  selectEdges(cw, nonEdges)
  deleteSelectedEdges(cw)
  setEdgeLineWidthDirect(cw, edges, edgeTable$width[edgeRows])
})

# Load new cytoscape session

# Term set overlap info
load("resources/termSetOverlap.RData")
overlapNumbers <- data.frame(sapply(termSetOverlap, function(r){sapply(r, length)}))
overlapNumbers <- new.order(pqPairs, overlapNumbers)
#max and min across regions
overlapNumbers[overlapNumbers < 10] <- NA
max(overlapNumbers, na.rm = T)
min(apply(overlapNumbers, 2, function(x) min(x, na.rm = TRUE)));
# Plot overlap
interaction_list2 <- apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  overlap <- as.numeric(overlapNumbers[, structure])
  width <- sapply(overlap, log1p)
  nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], color = avgExprColor[, structure])
  # nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], color = pQcolors)
  edgeTable <- data.frame(fromNode = pqPairs[ , 1], toNode = pqPairs[ , 2], overlap = overlap, 
                          edgeType="overlap", width = width)
  g <- cyPlot(nodeTable, edgeTable)
  windowName <- paste("overlapTermSets_",structName, sep = "")
  cw <- CytoscapeWindow(windowName, graph = g, overwrite = TRUE)
  displayGraph (cw)
  cyt.visuals(cw)
  setNodeColorDirect(cw,  nodeTable[ , "nodeName"],  nodeTable[ , "color"])
  
  edgeRows <- which(edgeTable$overlap >= 10)
  edges <- apply(edgeTable[edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (overlap) ")})
  nonEdges <- apply(edgeTable[-edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (overlap) ")})
  selectEdges(cw, nonEdges)
  deleteSelectedEdges(cw)
  setEdgeLineWidthDirect(cw, edges, edgeTable$width[edgeRows])
})
