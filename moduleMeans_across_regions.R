#Networks of module means in different regions to Cytoscape.
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

load("resources/polyQ.RData")
ontology <- read.csv("ABA_human_processed/Ontology_edited.csv")

###Cytoscape circular plots###
id <- ontology[ontology$name %in% "parietal lobe", ][ , c(1:3)]
load(paste("regional_coexpression/", gsub(" ", "_", id[3]), "/moduleMeans_", id[2], ".RData", sep = ""))
cyt <- exportNetworkToCytoscape(moduleMeans, threshold = 0.0, nodeAttr = pQcolors, 
          edgeFile = paste("regional_coexpression/", gsub(" ", "_", id[3]), "/moduleMeansEdges_", structName,".txt", sep = ""))

#########################################
#Variance of edges
regionLs <- split(structureIDs, seq(nrow(structureIDs)))
names(regionLs) <- gsub(" ", "_", structureIDs$name)
mm_list <- lapply(regionLs, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/moduleMeans_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  mm <- moduleMeans
  detach(2)
  mm
})

#Table with regional variance based on corr. values
pQpairs <- t(combn(polyQgenes, 2))
rowpairs <- apply(pQpairs, 1, function(x){paste( x[1], "-", x[2], sep = "")})
rsVar <- apply(simplify2array(mm_list), 1:2, var)
cyt <- exportNetworkToCytoscape(rsVar, edgeFile = "regional_variance/moduleMeansEdges_variance3.txt", threshold = 0.0)
regional_variance <- apply(pQpairs, 1, function(x){rsVar[x[1], x[2]]})
regions <- apply(pQpairs, 1, function(x){
  gene1 <- x[1]
  gene2 <- x[2]
  corVal <- sapply(mm_list, function(x){x[gene1, gene2]})
})
table1 <- cbind(regional_variance, t(regions))
rownames(table1) <- rowpairs
remove(regional_variance)
remove(regions)
table1 <- table1[order(table1[, "regional_variance"], decreasing = TRUE), ]

#Table with regional variance based on absolute corr. values
mm_list_abs <- lapply(mm_list, abs)
rsVar_abs <- apply(simplify2array(mm_list_abs), 1:2, var)
# cyt <- exportNetworkToCytoscape(rsVar_abs, edgeFile = "moduleMeansEdges_variance_abs.txt", threshold = 0.0)
regional_variance_abs <- apply(pQpairs, 1, function(x){rsVar_abs[x[1], x[2]]})
regions_abs <- apply(pQpairs, 1, function(x){
  gene1 <- x[1]
  gene2 <- x[2]
  corVal <- sapply(mm_list_abs, function(x){x[gene1, gene2]})
})
table2 <- cbind(regional_variance_abs, t(regions_abs))
rownames(table2) <- rowpairs
remove(regional_variance_abs)
remove(regions_abs)
table2 <- table2[order(table2[, "regional_variance_abs"], decreasing = TRUE), ]

#cat(paste("# Absolute correlation between polyQ genes across different regions.\n#", Sys.time(), "\n#gene_pair\t"), file = "table_correlations.txt")
#write.table(pQpairs, file = "table_correlations_abs.txt", sep = "\t", append = TRUE, quote = FALSE)
pdf(file = "table_correlations.pdf", 12, 16)
par(mar = c(6, 10, 15, 4));
labeledHeatmap(cbind(rep(0, 36), table1[, -1]), xLabels = colnames(table1), yLabels = rownames(table1), xLabelsPosition = "top", setStdMargins = FALSE,
               xLabelsAdj = 0, zlim = c(-1,1), colors = blueWhiteRed(200), textMatrix = round(table1, digits = 2),
               main = "Correlations between polyQ genes across different regions")
labeledHeatmap(cbind(rep(0, 36), table2[, -1]), xLabels = colnames(table2), yLabels = rownames(table2), xLabelsPosition = "top", setStdMargins = FALSE, 
               xLabelsAdj = 0, zlim = c(0,1), colors = blueWhiteRed(200)[100:200], textMatrix = round(table2, digits = 2), 
               main = "Absolute correlations between polyQ genes across different regions")
dev.off()
