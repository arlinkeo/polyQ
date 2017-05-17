# Plot co-expressed genes in HD-associated region
setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/polyQgenes")
library("rgexf")
options(stringsAsFactors = FALSE)

load("polyQ.RData")
load("regional_coexpression/HD_region/meanCor_HDnetworkBD.RData")
load("genesets_threshold050.RData")
probeInfo <- read.csv("/tudelft.net/staff-bulk/ewi/insy/DBL/sjoerdhuisman/ABA_human_brain_probegene/probe_info_2014-11-11.csv")
probeInfo[grep("TMEM189-UBE2V1", probeInfo$gene_name), 4] <- "TMEM189-UBE2V1" # change its gene symbol TMEM189 to TMEM189-UBE2V1, otherwise there are non-unique gene symbols
entrezId2Name <- function (x) { row <- match(x, probeInfo$entrez_id); probeInfo[row, 4]} #Input is single element

geneSets <- regionLs$HD_region
names(pQEntrezIDs) <- pQEntrezIDs
names(pQcolors) <- pQEntrezIDs
nodes <- unique(do.call(c, geneSets))
selection <- meanCor[nodes, nodes]
nodes <- data.frame(id = nodes, label = entrezId2Name(nodes))

membership <- as.data.frame(t(sapply(nodes$id, function(g){
  sapply(geneSets, function(set){
    g %in% set
  })
})))

nodeCol <- as.data.frame(t(sapply(nodes$id, function(g){
  ms <- unlist(membership[g, ])
  if (g %in% pQEntrezIDs) {
    col2rgb(pQcolors[g])
  }
  else if (length(ms) > 1){
    pqs <-  names(ms[which(ms)])
    colors <- pQcolors[pqs]
    mix.colors <- function(x){
      mixRgb <- apply(col2rgb(x), 1, mean)
      round(mixRgb)
    }
    mix.colors(colors)
  }
  else{#==1
    pq <-  names(ms[which(ms)])
    col2rgb(pQcolors[pq])
  }
})))
nodeCol$alpha <- rep(1, nrow(nodeCol))

colnames(membership) <- entrezId2Name(colnames(membership))
membership <- as.data.frame(apply(membership, c(1,2), as.numeric))

selection[lower.tri(selection)] <- 0
edges <- as.data.frame(as.table(selection))
colnames(edges) <- c("source", "target", "coexpr")
edges <- edges[which(edges$coexpr >0.5), ]

nodeSize <- data.frame(size = sapply(nodes$id, function(x){if (x %in% pQEntrezIDs) 50 else 20}))
nodeSize <- sapply(nodes$id, function(x){if (x %in% pQEntrezIDs) 50 else 20})
nodeViz <- list(size = nodeSize, color = nodeCol)

polyQ_gene <- as.integer(nodes$id %in% pQEntrezIDs)

write.gexf(nodes = nodes, edges = edges[, c("source", "target")], edgesWeight = edges[, c("coexpr")],
           edgesAtt = as.data.frame(edges[, c("coexpr")]), nodesVizAtt = nodeViz, nodesAtt = cbind(polyQ_gene, membership),
           defaultedgetype = "undirected", output = "HDregion.gexf")