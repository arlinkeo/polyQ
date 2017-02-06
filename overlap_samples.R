# Number samples of atomic structures that fall within the HD associated regions (Coppen2016)

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/BrainExpr.RData")
ontology <- read.csv("ABA_human_processed/Ontology_edited.csv")
load("resources/polyQ.RData")
structureIDs[, 3] <- sapply(structureIDs[, 3], function(id){gsub(" ", "_", id)})
rownames(structureIDs) <- structureIDs$name
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element

#Select HD region-specific samples
sampleIDs_HD <- lapply(donorNames, function(d){
  networkInfo <- read.csv(paste("HD_masks_Coppen2016/samples_in_networks_", d, ".txt", sep = ""), header = TRUE, sep = "\t", comment.char = "#")
  networkB <- networkInfo[, "network_B"]
  networkD <- networkInfo[, "network_D"]
  bitwOr(networkB, networkD)
})

#Select anatomic region-specific samples
sampleIDs <- apply(structureIDs, 1, function(id){
  print(id)
  structName <- id[2]
  ontologyRows <- grep(id[1], ontology$structure_id_path)
  selectIds <- as.character(ontology$id[ontologyRows])
  lapply(donorNames, function(d){
    expr <- brainExpr[[d]]
    ids <- intersect(selectIds, colnames(expr))
    cols <- colnames(expr) %in% ids
    print(paste("Samples: ", length(which(cols))))
    as.numeric(cols)
  })
})

sampleIDs <- c(sampleIDs, HD_region = list(sampleIDs_HD))

# Number of samples in atomic structures and HD regions
# overlap <- sapply(sampleIDs, function(s){
#   sapply(donorNames, function(d){
#     samples <- bitwAnd(s[[d]], sampleIDs_HD[[d]])
#     #paste(length(which(samples == 1)), "/", length(which(s[[d]] == 1)), sep = "")
#     #length(which(samples == 1)) / length(which(s[[d]] == 1)) * 100
#     #length(which(samples == 1)) / length(which(sampleIDs_HD[[d]] == 1)) * 100
#     length(which(samples == 1))
#   })
# })

# # Number of samples in structure
# sampleSize <- sapply(sampleIDs, function(s){
#   length(which(unlist(s) == 1))
# })

# Average expression of a polyQ gene in a structure across donors and samples.
avgExpr <- sapply(sampleIDs, function(s){
  res <- sapply(donorNames, function(d){
    expr <- brainExpr[[d]]
    expr2 <- expr[pQEntrezIDs, as.logical(s[[d]])]
    apply(expr2, 1, mean) # Avg across region-specific samples per PQ per donor
  })
  apply(res, 1, mean) # Avg across region-specific samples and donors per PQ
})
rownames(avgExpr) <- sapply(rownames(avgExpr), entrezId2Name)
pal <- colorRampPalette(c('darkgreen', 'yellow', 'red'))
avgExprColor <- matrix(pal(100)[as.numeric(cut(avgExpr, breaks = 100))], 9, 9)
rownames(avgExprColor) <- rownames(avgExpr)
colnames(avgExprColor) <- colnames(avgExpr)

load("resources/pqNeighbors.RData")

# Plot pqNeighbors in Cytoscape for each structure
pqPairs <- t(combn(polyQgenes, 2))
apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  mat <- pqNeighbors[[structure]]
  
  nodeTable  <- data.frame(nodeName = polyQgenes, expr = avgExpr[, structure], exprColor = avgExprColor[, structure])
  edgeTable <- data.frame(fromNode = pqPairs[ , 1], toNode = pqPairs[ , 2], edgeType="interaction",
                          interaction = apply(pqPairs, 1, function(x){mat[x[1], x[2]]}))
  g <- cyPlot(nodeTable, edgeTable)
  windowName <- paste("pqNeighbors_",structName, sep = "")
  cw <- CytoscapeWindow(windowName, graph = g, overwrite = TRUE)
  displayGraph (cw)
  layoutNetwork (cw, layout.name = "circular")
  setDefaultNodeShape (cw, "Ellipse")
  setDefaultNodeFontSize (cw, 8)
  setNodeColorDirect(cw,  nodeTable[ , "nodeName"],  nodeTable[ , "exprColor"])
  #lockNodeDimensions (cw, TRUE)
  setNodeHeightDirect(cw, nodeTable[ , "nodeName"], 30)
  setNodeWidthDirect(cw, nodeTable[ , "nodeName"], 30)
  #showGraphicsDetails(cw, TRUE)
  rowCol <- which(mat == 0, arr.ind = TRUE)
  nonEdges <- apply(rowCol, 1, function(x){
    row <- x[1]
    col <- x[2]
    paste(rownames(mat)[row], colnames(mat)[col], sep = " (interaction) ")
  })
  selectEdges(cw, nonEdges)
  deleteSelectedEdges(cw)
  svgName <- paste("regional_coexpression/", structure, "/", windowName, ".svg", sep = "")
  saveImage(cw, svgName, "svg")
  # cyt <- exportNetworkToCytoscape(pqNeighbors[[structure]], threshold = -1, nodeAttr = attr, 
  #                                 edgeFile = paste("regional_coexpression/", structure, "/pqNeighborsEdges_",structName,".txt", sep = ""), 
  #                               nodeFile = paste("regional_coexpression/", structure, "/pqExprNodes_", structName,".txt", sep = ""))
})

attr <- cbind(avgExpr[, "HD_region"], avgExprColor[, "HD_region"])
colnames(attr) <- c("expr", "exprColor")
cyt <- exportNetworkToCytoscape(pqNeighbors[["HD_region"]], threshold = -1, nodeAttr = attr, 
                                edgeFile = "HD_masks_Coppen2016/pqNeighborsEdges_HDregion.txt", 
                                nodeFile = "HD_masks_Coppen2016/pqExprNodes_HDregion.txt")
