# Plot co-expressed genes in HD-associated region
# setwd("/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/polyQgenes")
library("rgexf")

# load("polyQ.RData")
# load("regional_coexpression/HD_region/meanCor_HDnetworkBD.RData")
# load("genesets_threshold050.RData")
# probeInfo <- read.csv("/tudelft.net/staff-bulk/ewi/insy/DBL/sjoerdhuisman/ABA_human_brain_probegene/probe_info_2014-11-11.csv")
source("PolyQ_scripts/baseScript.R")
# load("regional_coexpression/HD_region/meanCor_HDnetworkBD.RData")
load("resources/genesets_threshold050.RData")
probeInfo[grep("TMEM189-UBE2V1", probeInfo$gene_name), 4] <- "TMEM189-UBE2V1" # change its gene symbol TMEM189 to TMEM189-UBE2V1, otherwise there are non-unique gene symbols

names(pQEntrezIDs) <- pQEntrezIDs
names(pQcolors) <- pQEntrezIDs

geneSets <- lapply(pQEntrezIDs, function(p) {
  v <-meanCor[p, ]
  top <- head(-sort(-v), sample(1:5, 1))#sample(10:20, 1))
  c(p, names(top[top > 0.5]))
})
# geneSets <- regionLs$HD_region
nodes <- unique(do.call(c, geneSets))

#random neighbours in network
nodes2 <- sample(setdiff(rownames(meanCor), nodes), length(nodes)/2) # nodes without pq's
nodes2 <- nodes2[sapply(nodes2, function(n) {
  v <-meanCor[n, nodes]
  if (sum(v > 0.5)) TRUE
  else FALSE
})]
nodes2 <- unique(do.call(c, as.list(nodes2)))
nodes <- unique(c(nodes, nodes2))

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
  else if (sum(ms) > 1){
    pqs <-  names(ms[which(ms)])
    colors <- pQcolors[pqs]
    mix.colors <- function(x){
      mixRgb <- apply(col2rgb(x), 1, mean)
      round(mixRgb)
    }
    mix.colors(colors)
  }
  else if (sum(ms) == 1) {#==1
    pq <-  names(ms[which(ms)])
    col2rgb(pQcolors[pq])
  }
  else{#==0
    col2rgb("#C0C0C0")
  }
})))
nodeCol$alpha <- rep(1, nrow(nodeCol))

colnames(membership) <- entrezId2Name(colnames(membership))
membership <- as.data.frame(apply(membership, c(1,2), as.numeric))

selection[lower.tri(selection)] <- 0
edges <- as.data.frame(as.table(selection))
colnames(edges) <- c("source", "target", "coexpr")
edges <- edges[which(edges$coexpr >0.5), ]

edgeCol <- t(apply(edges, 1, function(e){
  g1 <- e[[1]]
  g2 <- e[[2]]
  ms1 <- unlist(membership[g1, ])
  ms2 <- unlist(membership[g2, ])
  if (all(ms1 == ms2)) {
    unlist(nodeCol[g1, ])
  }
  else{#==0
    c(col2rgb("#C0C0C0"), 0.5)
  }
}))
edgeViz <- list(color = edgeCol)

nodeSize <- data.frame(size = sapply(nodes$id, function(x){if (x %in% pQEntrezIDs) 50 else 20}))
nodeSize <- sapply(nodes$id, function(x){if (x %in% pQEntrezIDs) 50 else 20})
nodeViz <- list(size = nodeSize, color = nodeCol)

polyQ_gene <- as.integer(nodes$id %in% pQEntrezIDs)

write.gexf(nodes = nodes, edges = edges[, c("source", "target")], edgesWeight = edges[, c("coexpr")], edgesVizAtt = edgeViz, 
           edgesAtt = as.data.frame(edges[, c("coexpr")]), nodesVizAtt = nodeViz, nodesAtt = cbind(polyQ_gene, membership),
           defaultedgetype = "undirected", output = "Images/HDregion_network/HDregion.gexf")