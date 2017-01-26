#Search for polyQ genes in other polyQ gene sets

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
allowWGCNAThreads(nThreads = 32)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
name2entrezId <- function (x) { row <- which(probeInfo$gene_symbol == x); probeInfo[row, 6]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

#Load gene sets (EntrezIDs)
load("resources/genesets_threshold050.RData")
load("resources/genesets_threshold050_HDregion.RData")
regionLs <- c(regionLs, HD_region = list(selection))

genepairs <- t(combn(polyQgenes, 2))
rownames(genepairs) <- apply(genepairs, 1, function(x) {paste(x[1], "-", x[2], sep = "") })


reduce2pq <- function(r){lapply(r, function(x){intersect(x, pQEntrezIDs)})}
pqInteractors <- lapply(regionLs, reduce2pq)
pqInteractors <- lapply(pqInteractors, function(r){r[which(lengths(r) != 1)]}) # remove vectors of size 1
pqInteractors <- pqInteractors[which(lengths(pqInteractors) != 0)] # remove empty lists
pqInteractors <- lapply(pqInteractors, function(r){unique(unlist(r))})
pqInteractors <- lapply(pqInteractors, function(r){lapply(r, entrezId2Name)})

lapply(pqInteractors, function(x){unlist(lapply(names(x), entrezId2Name))})
pqInteractors <- 