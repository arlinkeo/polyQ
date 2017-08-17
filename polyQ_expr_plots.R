# Distribution of gene expression
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
#library(WGCNA)
options(stringsAsFactors = FALSE)

load("resources/sampleIDs.RData")
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("brain", "cerebellum"), ]
structureIDs <- rbind(HD_region = c(NA, "HDregion", "HD_region"), structureIDs)
sampleIDs <- sampleIDs[!names(sampleIDs) %in% c("brain", "cerebellum")]
load("resources/brainExpr.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]}
ontology <- read.csv("ABA_human_processed/Ontology_edited.csv")
rownames(ontology) <- ontology$id

pq_expr <- lapply(brainExpr, function(mat){mat[pQEntrezIDs, ]})

###############################################################################
# #Reduce samples to those that occur in all donors, and convert row- and colnames
# all_samples <- lapply(donorList, colnames)
# unique_samples <- unique(unlist(all_samples))
# alldonor_samples <- unique_samples[sapply(unique_samples, function(x){if (FALSE %in% sapply(all_samples, function(y){x %in% y})) FALSE else TRUE})]
# alldonor_samples <- names(sort(sapply(alldonor_samples, function(x){ontology[x, 'graph_order']})))
# donorList <- lapply(donorList, function(x){x[ , alldonor_samples]})
# donorList <- lapply(donorList, function(x){
#   colnames(x) <- sapply(colnames(x), function(y){ontology[y, 'acronym']})
#   x
# })
# 
# #Mean and std. dev. for each polyQ gene in all samples across donors
# mean_expr <- apply(simplify2array(donorList), 1:2, mean)
# sd_expr <- apply(simplify2array(donorList), 1:2, sd)
# #var_expr <- apply(simplify2array(donorList), 1:2, var)
# 
# #Plot histogram
# colors <- sapply(alldonor_samples, function(x){paste("#", ontology[x, 'color_hex_triplet'], sep = "")})
# pdf(file = "polyQ_expr_plots.pdf", 24)
# for (pq in polyQgenes){
#   lapply(names(donorList), function(d){
#     pq_expr <- unlist(donorList[[d]][pq, ])
#     barplot(pq_expr, main = paste(pq, " expression in donor ", d, sep = ""), las = 2, col = colors, ylim = c(0,10))
#   })
#   barplot(mean_expr[pq, ], main = paste("Average ", pq, " expression across donors ", sep = ""), las = 2, col = colors, ylim = c(0,10))
#   barplot(sd_expr[pq, ], main = paste("Standard deviation of ", pq, " expression across donors ", sep = ""), las = 2, col = colors, ylim = c(0,2.5))
# }
# dev.off()

############################################################################
### Plot per donor all samples ###

pdf(file = "polyQ_expr_plots5.pdf", 60, 40)

lapply(donorNames, function(d){
  expr <- as.matrix(pq_expr[[d]])
  colOrder <- names(sort(sapply(colnames(expr), function(x) {ontology[x, 'graph_order']})))
  expr <- expr[ , colOrder]
  colors <- sapply(colnames(expr), function(x){
    clr <- ontology[x, 'color_hex_triplet']
    if (nchar(clr) == 5) {paste("#0", clr, sep = "")}
    else {paste("#", clr, sep = "")}
  })
  colnames(expr) <- sapply(colnames(expr), function(x){ontology[x, 'acronym']})
  par(mfrow = c(9, 1), oma = c(0, 0, 15, 0), lwd = 0.01, mai = c(1, 0, 1, 0));
  sapply(rownames(expr), function(x){
    barplot(expr[x, ], main = bquote(italic(.(entrezId2Name(x)))), cex.main = 10, col = colors, las = 2, cex.names = 0.5, 
            ylim = c(0, 10.5), yaxt = 'n', border = "gray")
    axis(2, pos = 0, cex.axis = 3)
  })
  title(paste("Donor ", d, sep = ""), outer = TRUE, cex.main = 8)
})

dev.off()

############################################################################
### Plot per donor all samples from all regions
# Plot
pdf(file = "polyQ_expr_plots3.pdf", 80, 30)#60, 40) #80,30

lapply(donorNames, function(d){
  expr <- as.matrix(pq_expr[[d]])
  ids <- lapply(sampleIDs[-1], function(r){v <- r[[d]]})
  ids <- as.logical(Reduce('|', ids))
  expr <- expr[ , ids]
  colOrder <- names(sort(sapply(colnames(expr), function(x) {ontology[x, 'graph_order']})))
  expr <- expr[ , colOrder]
  colors <- sapply(colnames(expr), function(x){
    clr <- ontology[x, 'color_hex_triplet']
    if (nchar(clr) == 5) {paste("#0", clr, sep = "")}
    else {paste("#", clr, sep = "")}
  })
  colnames(expr) <- sapply(colnames(expr), function(x){ontology[x, 'acronym']})
  par(mfrow = c(9, 1), oma = c(0, 0, 15, 0), lwd = 0.01);
  sapply(rownames(expr), function(x){
    barplot(expr[x, ], main = bquote(italic(.(entrezId2Name(x)))), cex.main = 8, col = colors, las = 2, cex.names = 0.5, 
            ylim = c(0, 10.5), yaxt = 'n')
    axis(2, pos = 0, cex.axis = 3)
  })
  title(paste("Donor ", d, sep = ""), outer = TRUE, cex.main = 8)
})

dev.off()