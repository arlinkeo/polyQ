# Plot data in Cytoscape

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(RCy3)
options(stringsAsFactors = FALSE)

load("resources/polyQ.RData")
structureIDs[, 3] <- sapply(structureIDs[, 3], function(id){gsub(" ", "_", id)})
structureIDs <- rbind(structureIDs, c(NA, "HDregion", "HD_region"))
rownames(structureIDs) <- structureIDs$name
pqPairs <- t(combn(polyQgenes, 2))
rownames(pqPairs) <- apply(pqPairs, 1, function(x){paste(x[1], "-", x[2], sep = "")})

# Node info
load("resources/avgExpr.RData")
load("resources/avgExprColor.RData")

linMap <- function(x){
  if (x>0.5) (x - 0.5) * 2 * 10 + 2
  else 0
}
# Cytoscape default visuals
cyt.visuals <- function(cw) {
  layoutNetwork (cw, layout.name = "circular")
  setDefaultNodeShape (cw, "Ellipse")
  setDefaultNodeFontSize (cw, 8)
  setNodeHeightDirect(cw, polyQgenes, 30)
  setNodeWidthDirect(cw, polyQgenes, 30)
  #lockNodeDimensions (cw, TRUE) # Not implemented yet in RCy3
  #setDefaultNodeSize(cw, 30)
  #showGraphicsDetails(cw, TRUE) # Not implemented yet in RCy3
}

# Load info about direct interacters (single_corr.R)
load("resources/sc_list.RData")
# Plot pqNeighbors in Cytoscape for each structure
apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  mat <- sc_list[[structure]] # adjacency matrix
  width <- apply(mat, c(1,2), linMap)
  
  nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], exprColor = avgExprColor[, structure])
  edgeTable <- data.frame(fromNode = pqPairs[ , 1], toNode = pqPairs[ , 2], edgeType="interaction",
                          interaction = apply(pqPairs, 1, function(x){mat[x[1], x[2]]}), 
                          width = apply(pqPairs, 1, function(x){width[x[1], x[2]]}))
  g <- cyPlot(nodeTable, edgeTable)
  windowName <- paste("pqNeighbors_",structName, sep = "")
  cw <- CytoscapeWindow(windowName, graph = g, overwrite = TRUE)
  displayGraph (cw)
  cyt.visuals(cw)
  setNodeColorDirect(cw,  nodeTable[ , "nodeName"],  nodeTable[ , "exprColor"])
  edgeRows <- which(edgeTable$interaction > 0.5)
  edges <- apply(edgeTable[edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (interaction) ")})
  nonEdges <- apply(edgeTable[-edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (interaction) ")})
  selectEdges(cw, nonEdges)
  deleteSelectedEdges(cw)
  setEdgeLineWidthDirect(cw, edges, edgeTable$width[edgeRows])
  
  # svgName <- paste("Images/pqNeighbors2_", structure, "/", windowName, ".svg", sep = "")
  # saveImage(cw, svgName, "svg")
})

# Load new session

# Gene set overlap info
load("resources/geneSetOverlap.RData")
overlapNumbers <- sapply(geneSetOverlap, function(r){sapply(r, length)})
overlapNumbers <- overlapNumbers[rownames(pqPairs), ]
load("resources/geneSetOverlapSignif.RData")
# Plot for overlap of gene sets
interaction_list <- apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  print(structure)
  overlap <- overlapNumbers[, structure]
  width <- sapply(overlap, log1p)
  signif <- geneSetOverlapSignif[, structure] # adjacency matrix
  
  nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], exprColor = avgExprColor[, structure])
  edgeTable <- data.frame(fromNode = pqPairs[ , 1], toNode = pqPairs[ , 2], edgeType="overlap",
                          overlap = overlap, width = width, signif = signif)
  g <- cyPlot(nodeTable, edgeTable)
  windowName <- paste("overlapGeneSets_",structName, sep = "")
  cw <- CytoscapeWindow(windowName, graph = g, overwrite = TRUE)
  displayGraph (cw)
  cyt.visuals(cw)
  setNodeColorDirect(cw,  nodeTable[ , "nodeName"],  nodeTable[ , "exprColor"])
  
  edgeRows <- which(edgeTable$signif < 0.05)
  edges <- apply(edgeTable[edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (overlap) ")})
  nonEdges <- apply(edgeTable[-edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (overlap) ")})
  selectEdges(cw, nonEdges)
  deleteSelectedEdges(cw)
  setEdgeLineWidthDirect(cw, edges, edgeTable$width[edgeRows])
  
  # svgName <- paste("regional_coexpression/", structure, "/", windowName, ".svg", sep = "")
  # saveImage(cw, svgName, "svg")
})

# Term set overlap info
load("resources/termSetOverlap.RData")
#load("resources/termSetOverlapSignif.RData")
# Plot overlap
interaction_list2 <- apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  mat <- termSetOverlap[[structure]] # adjacency matrix
  width <- apply(mat, c(1,2), log1p)
  #matSignif <- termSetOverlapSignif[[structure]] # adjacency matrix
  
  nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], exprColor = avgExprColor[, structure])
  edgeTable <- data.frame(fromNode = pqPairs[ , 1], toNode = pqPairs[ , 2], edgeType="overlap",
                          overlap = apply(pqPairs, 1, function(x){mat[x[1], x[2]]}), 
                          width = apply(pqPairs, 1, function(x){width[x[1], x[2]]})) 
                          #overlapSignif = apply(pqPairs, 1, function(x){matSignif[x[1], x[2]]}))
  g <- cyPlot(nodeTable, edgeTable)
  windowName <- paste("overlapTermSets_",structName, sep = "")
  cw <- CytoscapeWindow(windowName, graph = g, overwrite = TRUE)
  displayGraph (cw)
  cyt.visuals(cw)
  setNodeColorDirect(cw,  nodeTable[ , "nodeName"],  nodeTable[ , "exprColor"])
  
  edgeRows <- which(edgeTable$overlap > 10)
  edges <- apply(edgeTable[edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (overlap) ")})
  nonEdges <- apply(edgeTable[-edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (overlap) ")})
  selectEdges(cw, nonEdges)
  deleteSelectedEdges(cw)
  setEdgeLineWidthDirect(cw, edges, edgeTable$width[edgeRows])

  # svgName <- paste("regional_coexpression/", structure, "/", windowName, ".svg", sep = "")
  # saveImage(cw, svgName, "svg")
})

#Plot mean coexpression between top 25 gene sets
apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  file <- paste("regional_coexpression/", structure, "/moduleMeans_", structName, ".RData", sep = "")
  attach(file)
  mat <- moduleMeans
  detach(2)
  
  # pal <- colorRampPalette(c('red', 'white', 'blue'))
  # mmColor <- matrix(pal(200)[as.numeric(cut(mat, breaks = 100))], 9, 9)
  # rownames(mmColor) <- rownames(mat)
  # colnames(mmColor) <- colnames(mat)
  width <- apply(mat, c(1,2), linMap)
  
  nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], exprColor = avgExprColor[, structure])
  edgeTable <- data.frame(fromNode = pqPairs[ , 1], toNode = pqPairs[ , 2], edgeType="coexpr",
                          coexpr = apply(pqPairs, 1, function(x){mat[x[1], x[2]]}), 
                          # color = apply(pqPairs, 1, function(x){mmColor[x[1], x[2]]}), 
                          width = apply(pqPairs, 1, function(x){width[x[1], x[2]]}))
  g <- cyPlot(nodeTable, edgeTable)
  windowName <- paste("meanCoexpr_",structName, sep = "")
  cw <- CytoscapeWindow(windowName, graph = g, overwrite = TRUE)
  displayGraph (cw)
  cyt.visuals(cw)
  setNodeColorDirect(cw,  nodeTable[ , "nodeName"],  nodeTable[ , "exprColor"])
  
  edgeRows <- which(edgeTable$coexpr > 0.5)
  edges <- apply(edgeTable[edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (coexpr) ")})
  nonEdges <- apply(edgeTable[-edgeRows, ], 1, function(e){paste(e[1], e[2], sep = " (coexpr) ")})
  selectEdges(cw, nonEdges)
  deleteSelectedEdges(cw)
  setEdgeLineWidthDirect(cw, edges, edgeTable$width[edgeRows])
  #setEdgeColorDirect(cw, apply(edgeTable, 1, function(x){paste(x[1], x[2], sep = " (coexpr) ")}), edgeTable$color)

  # svgName <- paste("regional_coexpression/", structure, "/", windowName, ".svg", sep = "")
  # saveImage(cw, svgName, "svg")
})