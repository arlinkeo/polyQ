# Plot data in Gephi
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library("rgexf")
options(stringsAsFactors = FALSE)

load("resources/polyQ.RData")
load("resources/avgExpr.RData")
load("resources/avgExprColor.RData")

structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellum"), ]
structureIDs <- rbind(c(NA, "HDregion", "HD_region"), structureIDs)
rownames(structureIDs) <- structureIDs$name

pqOrder <- data.frame(order = as.integer(c(c(1:3), 9, c(4:8))))
nodes <- data.frame(id = pQEntrezIDs, label = polyQgenes)
nodeSize <- data.frame(size = rep(50, length(pQcolors)))
rownames(nodeSize) <- polyQgenes

structures <- structureIDs$name
names(structures) <- structureIDs$name
nodeViz <- lapply(structures, function(s){
  nodeCol <- data.frame(t(col2rgb(avgExprColor[, s])), alpha = rep(1, length(polyQgenes)))
  list(size = nodeSize, color = nodeCol)
})

# Load info about direct interacters (single_corr.R)
load("resources/sc_list.RData")
# Plot pqNeighbors gexf-file in Gephi for each structure
apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  fName <- paste("Images/", "pqNeighbors_", structName, ".gexf", sep = "")
  mat <- sc_list[[structure]] # adjacency matrix
  rownames(mat) <- pQEntrezIDs
  colnames(mat) <- pQEntrezIDs
  mat[lower.tri(mat)] <- 0
  edges <- as.data.frame(as.table(mat))
  colnames(edges) <- c("source", "target", "coexpr")
  rows <- which(edges$coexpr >0.5)
  if (length(rows) != 0) edges <- edges[rows, ]

  write.gexf(nodes = nodes, edges = edges[, c("source", "target")], edgesWeight = edges[, c("coexpr")],
              nodesAtt = pqOrder, nodesVizAtt = nodeViz[[structure]], 
             defaultedgetype = "undirected", output = fName)
})

# Gene set overlap info
load("resources/geneSetOverlap.RData")
overlapNumbers <- sapply(geneSetOverlap, function(r){sapply(r, length)})
load("resources/geneSetOverlapSignif.RData")
# Plot for overlap of gene sets
apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  fName <- paste("Images/", "overlapGeneSets_", structName, ".gexf", sep = "")
  overlap <- overlapNumbers[, structure]
  width <- sapply(overlap, log1p)
  signif <- geneSetOverlapSignif[, structure]
  pairs <- t(sapply(names(overlap), function(pair){v <- unlist(strsplit(pair, split = "-")); nodes$id[nodes$label %in% v]}))
  edges <- data.frame(source = pairs[ ,1], target = pairs[ ,2], overlap, width)
  rows <- which(signif < 0.05)
  if (length(rows) != 0) edges <- edges[rows, ]

  write.gexf(nodes = nodes, edges = edges[, c("source", "target")], edgesWeight = edges[, c("width")],
             nodesAtt = pqOrder, nodesVizAtt = nodeViz[[structure]], edgesAtt = data.frame(overlap = edges[, c("overlap")]),
             defaultedgetype = "undirected", output = fName)
})

# Term set overlap info
load("resources/termSetOverlap.RData")
# Plot overlap
apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  fName <- paste("Images/", "overlapTermSets_", structName, ".gexf", sep = "")
  overlap <- sapply(termSetOverlap[[structure]], length)
  width <- sapply(overlap, log1p)
  pairs <- t(sapply(names(overlap), function(pair){v <- unlist(strsplit(pair, split = "-")); nodes$id[nodes$label %in% v]}))
  edges <- data.frame(source = pairs[ ,1], target = pairs[ ,2], overlap, width)
  rows <- which(edges$overlap >= 10)
  if (length(rows) != 0) edges <- edges[rows, ]

  write.gexf(nodes = nodes, edges = edges[, c("source", "target")], edgesWeight = edges[, c("overlap")],
             nodesAtt = pqOrder, nodesVizAtt = nodeViz[[structure]], edgesAtt = data.frame(overlap = edges[, c("overlap")]),
             defaultedgetype = "undirected", output = fName)
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