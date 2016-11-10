# Check for overlap in GO terms of gene sets per brain region

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
#library(WGCNA)
library("RDAVIDWebService")
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellar nuclei","basal forebrain","globus pallidus"), ] # remove structures from list
structureIDs[, 3] <- sapply(structureIDs[, 3], function(x){gsub(" ", "_", x)})
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}
region.acronym <- function(x) {structureIDs[structureIDs$name %in% x, ]$acronym}

genepairs <- t(combn(pQEntrezIDs, 2))
rownames(genepairs) <- apply(genepairs, 1, function(x){paste(entrezId2Name(x[1]), "-", entrezId2Name(x[2]), sep = "")})