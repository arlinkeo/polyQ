# Check for overlap in polyQ gene sets per region

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellar nuclei","basal forebrain","globus pallidus"), ] # remove structures from list
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

genepairs <- t(combn(polyQgenes, 2))
rownames(genepairs) <- apply(genepairs, 1, function(x){paste(x[1], "-", x[2], sep = "")})

# Function lists duplicates in one region
duplicates <- function(x){
  pQidx <- which(colnames(x) %in% pQEntrezIDs)
  idx_dups <- which(duplicated(colnames(x)))
  dups_id <- colnames(x)[idx_dups]
  bin_matrix <- matrix(0, length(dups_id), length(polyQgenes), dimnames = list(dups_id, polyQgenes))
  for (id in dups_id) { # for each duplicated gene
    col_idx <- which(colnames(x) == id)
    pqsets <- sapply(col_idx, function(y){ polyQgenes[tail(which(pQidx < y), n = 1)] })
    for (pq in pqsets) {
      bin_matrix[id, pq] <- 1
    }
  }
  rownames(bin_matrix) <- lapply(rownames(bin_matrix), entrezId2Name)
  bin_matrix
}

# Function to count number of overlapping genes between two polyQ sets
mk.table <- function(ll){
  sapply(ll, function(r){
    apply(genepairs, 1, function(x){
      pq1 <- r[ , x[1]]
      pq2 <- r[ , x[2]]
      sum(bitwAnd(pq1, pq2))
    })
  })
}



attach("resources/regionLs_threshold080.RData")
duplicates_080 <- lapply(regionLs, function(x){duplicates(x)})
detach(2)
attach("resources/regionLs_threshold070.RData")
duplicates_070 <- lapply(regionLs, function(x){duplicates(x)})
detach(2)
attach("resources/regionLs_threshold060.RData")
duplicates_060 <- lapply(regionLs, function(x){duplicates(x)})
detach(2)

### Count number of overlapping genes between two polyQ sets ###
table_080 <- mk.table(duplicates_080)
table_070 <- mk.table(duplicates_070)
table_060 <- mk.table(duplicates_060)

#Load asssociations info from literature and combine with overlap table
asssociations <- read.csv("datatype_interactions.txt", sep = "\t", row.names = 1, comment.char = "#")
asssociations <- asssociations[, ncol(asssociations), drop = FALSE]
asssociations <- asssociations[rownames(genepairs), ]
table_080 <- cbind(table_080, asssociations)
table_070 <- cbind(table_070, asssociations)
table_060 <- cbind(table_060, asssociations)

# Plot table
pdf(file = "overlap_genesets.pdf", 12, 16)
par(mar = c(6, 10, 15, 4));
labeledHeatmap(as.matrix((table_060 > 0) + 0), xLabels = colnames(table_060), yLabels = make.italic(rownames(table_060)), setStdMargins = FALSE,
               xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), textMatrix = table_060, 
               main = "Overlap between two polyQ gene sets with genes correlated >0.6")
labeledHeatmap(as.matrix((table_070 > 0) + 0), xLabels = colnames(table_070), yLabels = make.italic(rownames(table_070)), setStdMargins = FALSE,
               xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), textMatrix = table_070, 
               main = "Overlap between two polyQ gene sets with genes correlated >0.7")
labeledHeatmap(as.matrix((table_080 > 0) + 0), xLabels = colnames(table_080), yLabels = make.italic(rownames(table_080)), setStdMargins = FALSE,
               xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), textMatrix = table_080, 
               main = "Overlap between two polyQ gene sets with genes correlated >0.8")
dev.off()