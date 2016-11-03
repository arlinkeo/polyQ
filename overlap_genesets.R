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

load("resources/regionLs_threshold070.RData")

# Function lists duplicates in one region
duplicates <- function(x){
  pQidx <- which(colnames(x) %in% pQEntrezIDs)
  idx_dups <- which(duplicated(colnames(x)))
  dups_id <- colnames(x)[idx_dups]
  bin_matrix <- matrix(0, length(dups_id), length(polyQgenes), dimnames = list(dups_id, polyQgenes))
  for (id in dups_id) { # for each duplicated gene, output list of polyQ sets
    col_idx <- which(colnames(x) == id)
    pqsets <- sapply(col_idx, function(y){ polyQgenes[tail(which(pQidx < y), n = 1)] })
    for (pq in pqsets) {
      bin_matrix[id, pq] <- 1
    }
  }
  rownames(bin_matrix) <- lapply(rownames(bin_matrix), entrezId2Name)
  bin_matrix
}

duplicates_070 <- lapply(regionLs, function(x){duplicates(x)})

### Count number of overlapping genes between two polyQ sets ###
genepairs <- t(combn(polyQgenes, 2))
rownames(genepairs) <- apply(genepairs, 1, function(x){paste(x[1], "-", x[2], sep = "")})

table <- sapply(duplicates_070, function(r){
  apply(genepairs, 1, function(x){
    pq1 <- r[ , x[1]]
    pq2 <- r[ , x[2]]
    sum(bitwAnd(pq1, pq2))
  })
})

