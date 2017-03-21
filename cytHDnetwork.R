# Plot co-expressed genes in HD-associated region
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(RCy3)
options(stringsAsFactors = FALSE)

load("resources/polyQ.RData")
load("regional_coexpression/HD_region/meanCor_HDnetworkBD.RData")
load("resources/genesets_threshold050.RData")
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element

geneSets <- regionLs$HD_region
names(pQEntrezIDs) <- pQEntrezIDs
geneSets <- lapply(pQEntrezIDs, function(pq){
  set <- geneSets[[pq]]
  mat <- meanCor[pq, set[-1]]
  mat <- -sort(-mat)
  top50 <- head(mat, 50)
  c(set[1], names(top50))
})

interactors <- unique(do.call(c, geneSets))
selection <- meanCor[interactors, interactors]
nodeNames <- sapply(interactors, entrezId2Name)
rownames(selection) <- nodeNames
colnames(selection) <- nodeNames
names(geneSets) <- sapply(names(geneSets), entrezId2Name)
geneSets <- lapply(geneSets, function(set){sapply(set, entrezId2Name)})

mix.colors <- function(x){
  mixRgb <- apply(col2rgb(x), 1, mean)
  rgb(mixRgb[1], mixRgb[2], mixRgb[3], maxColorValue = 255)
}

nodeCol <- sapply(nodeNames, function(g){
  membership <- sapply(geneSets, function(set){
    g %in% set
  })
  if (length(membership) > 1){
    pqs <-  names(membership[which(membership)])
    colors <- pQcolors[pqs]
    mix.colors(colors)
  }
  else{#==1
    pq <-  names(membership[which(membership)])
    pQcolors[pq]
  }
})
names(nodeCol) <- nodeNames

labelSize <- sapply(nodeNames, function(x){if (x %in% polyQgenes) 16 else 8})
nodeSize <- sapply(nodeNames, function(x){if (x %in% polyQgenes) 80 else 30})
coexpr <- as.data.frame(as.table(selection))

#plot graph
nodeTable  <- data.frame(nodeName = nodeNames, nodeAttr =rep(1, length(nodeNames)), color = nodeCol, 
                         labelSize = labelSize, nodeSize = nodeSize)
edgeTable <- data.frame(fromNode = coexpr[ , 1], toNode = coexpr[ , 2], edgeType="coexpr",
                        coexpr = coexpr[ , 3])#, 
                        #width = apply(pqPairs, 1, function(x){width[x[1], x[2]]}))
edgeTable <- edgeTable[which(edgeTable$coexpr >0.5), ]
g <- cyPlot(nodeTable, edgeTable)
windowName <- "HDregion050"
cw <- CytoscapeWindow(windowName, graph = g, overwrite = TRUE)
displayGraph (cw)

layoutNetwork (cw, layout.name = "attribute-circle")
setDefaultNodeShape (cw, "Ellipse")
#setDefaultNodeFontSize (cw, 8)
setNodeFontSizeDirect(cw, nodeNames, labelSize)
setNodeSizeDirect(cw, nodeNames, nodeSize)
setNodeColorDirect(cw,  nodeTable[ , "nodeName"],  nodeTable[ , "color"])
#lockNodeDimensions (cw, TRUE) # Not implemented yet in RCy3
#setDefaultNodeSize(cw, 30)
#showGraphicsDetails(cw, TRUE) # Not implemented yet in RCy3

lapply(unique(nodeTable$color), function(c){ # c= unique(nodeTable$color)[6]
  nodes <- nodeNames[which(nodeTable$color %in% c)]
  selectNodes(cw, nodes)
  clearSelection(cw)
})