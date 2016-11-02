# Distribution of gene expression

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Load gene ID's
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]}
load("resources/polyQ.RData")
#probeInfo[probeInfo$entrez_id %in% pQEntrezIDs, ]
genes <- probeInfo[ , 6]
ontology <- read.csv("ABA_human_processed/Ontology_edited.csv")
rownames(ontology) <- ontology$id

#Read expression data
# donorList <- list("9861" = "9861",
#                   "10021" = "10021", 
#                   "12876" = "12876", 
#                   "14380" = "14380", 
#                   "15496" = "15496", 
#                   "15697" = "15697")
# donorList <- lapply(donorList, function(x){
#   gene_expr <- read.csv(paste("ABA_human_processed/gene_expr_normalized_microarray_donor", x, "_2014-11-11.csv", sep = ""), header = FALSE)
#   rownames(gene_expr) <- genes
#   sample_info <- read.csv(paste("ABA_human_processed/sample_info_normalized_microarray_donor", x, "_2014-11-11.csv", sep = ""))
#   colnames(gene_expr) <- sample_info$structure_id
#   gene_expr <- gene_expr[pQEntrezIDs, ]
#   rownames(gene_expr) <- sapply(rownames(gene_expr), entrezId2Name)
#   as.matrix(gene_expr)
# })
# remove(genes)

load("resources/polyQ_expr.RData")

###############################################################################
#Reduce samples to those that occur in all donors, and convert row- and colnames
all_samples <- lapply(donorList, colnames)
unique_samples <- unique(unlist(all_samples))
alldonor_samples <- unique_samples[sapply(unique_samples, function(x){if (FALSE %in% sapply(all_samples, function(y){x %in% y})) FALSE else TRUE})]
alldonor_samples <- names(sort(sapply(alldonor_samples, function(x){ontology[x, 'graph_order']})))
donorList <- lapply(donorList, function(x){x[ , alldonor_samples]})
donorList <- lapply(donorList, function(x){
  colnames(x) <- sapply(colnames(x), function(y){ontology[y, 'acronym']})
  x
})

#Mean and std. dev. for each polyQ gene in all samples across donors
mean_expr <- apply(simplify2array(donorList), 1:2, mean)
sd_expr <- apply(simplify2array(donorList), 1:2, sd)
#var_expr <- apply(simplify2array(donorList), 1:2, var)

#Plot histogram
colors <- sapply(alldonor_samples, function(x){paste("#", ontology[x, 'color_hex_triplet'], sep = "")})
pdf(file = "polyQ_expr_plots.pdf", 24)
for (pq in polyQgenes){
  lapply(names(donorList), function(d){
    pq_expr <- unlist(donorList[[d]][pq, ])
    barplot(pq_expr, main = paste(pq, " expression in donor ", d, sep = ""), las = 2, col = colors, ylim = c(0,10))
  })
  barplot(mean_expr[pq, ], main = paste("Average ", pq, " expression across donors ", sep = ""), las = 2, col = colors, ylim = c(0,10))
  barplot(sd_expr[pq, ], main = paste("Standard deviation of ", pq, " expression across donors ", sep = ""), las = 2, col = colors, ylim = c(0,2.5))
}
dev.off()

############################################################################
### Plot per donor all samples ###

pdf(file = "polyQ_expr_plots2.pdf", 60, 40) #80,30

lapply(names(donorList), function(d){
  expr <- donorList[[d]]
  expr <- expr[ , names(sort(sapply(colnames(expr), function(x) {ontology[x, 'graph_order']})))]
  colors <- sapply(colnames(expr), function(x){
    clr <- ontology[x, 'color_hex_triplet']
    if (nchar(clr) == 5) {paste("#0", clr, sep = "")}
    else {paste("#", clr, sep = "")}
  })
  colnames(expr) <- sapply(colnames(expr), function(x){ontology[x, 'acronym']})
  par(mfrow = c(9, 1), oma = c(0, 0, 15, 0), lwd = 0.01);
  sapply(rownames(expr), function(x){
    barplot(expr[x, ], main = bquote(italic(.(x))), cex.main = 5, col = colors, las = 2, cex.names = 0.5, 
            ylim = c(0, 10.5), yaxt = 'n')
    axis(2, pos = 0, cex.axis = 3)
  })
  title(paste("Donor ", d, sep = ""), outer = TRUE, cex.main = 8)
})

dev.off()

############################################################################
### Plot per donor all samples from 6 regions: cerebellar cortex, striatum, mesencephalon, pons, hypothalamus, and frontal lobe. ###
### whole brain is excluded ###
structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellar nuclei","basal forebrain","globus pallidus", "brain"), ] # remove structures from list

#Select all region-specific samples
sampleIds <- unlist(apply(structureIDs, 1, function(x){
  ontologyRows <- grep(x[1], ontology$structure_id_path)
  selectIds <- as.character(ontology$id[ontologyRows])
  print(paste(length(selectIds), "selected structure id's in", x[3]))
  selectIds
}))

# Plot
pdf(file = "polyQ_expr_plots3.pdf", 60, 40) #80,30

lapply(names(donorList), function(d){
  expr <- donorList[[d]]
  ids <- intersect(sampleIds, colnames(expr))
  cols <- colnames(expr) %in% ids
  expr <- expr[ , cols]
  print(paste("Samples: ", length(which(cols)), "in donor", d))
  expr <- expr[ , names(sort(sapply(colnames(expr), function(x) {ontology[x, 'graph_order']})))]
  colors <- sapply(colnames(expr), function(x){
    clr <- ontology[x, 'color_hex_triplet']
    if (nchar(clr) == 5) {paste("#0", clr, sep = "")}
    else {paste("#", clr, sep = "")}
  })
  colnames(expr) <- sapply(colnames(expr), function(x){ontology[x, 'acronym']})
  par(mfrow = c(9, 1), oma = c(0, 0, 15, 0), lwd = 0.01);
  sapply(rownames(expr), function(x){
    barplot(expr[x, ], main = bquote(italic(.(x))), cex.main = 5, col = colors, las = 2, cex.names = 0.5, 
            ylim = c(0, 10.5), yaxt = 'n')
    axis(2, pos = 0, cex.axis = 3)
  })
  title(paste("Donor ", d, sep = ""), outer = TRUE, cex.main = 8)
})

dev.off()