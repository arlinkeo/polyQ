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
load("resources/genesets_threshold070.RData")
region.acronym <- function(x) {structureIDs[structureIDs$name %in% x, ]$acronym}
# Open connection to DAVID
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl", url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

# Obtain a list of GO terms for each gene set
t <- 0.05 # EASE p-value threshold
lapply(names(regionLs), function(r){
  lapply(regionLs[[r]], function(pq){
    if (length(pq) > 1){
      pqname <- entrezId2Name(pq[1])
      listname <- paste(r, "_", pqname, sep = "")
      # Backgroundlist?
      result <- addList(david, pq, idType = "ENTREZ_GENE_ID", listName = listname, listType = "Gene")
      fname <- paste("regional_coexpression/", r, "/goterms070_", region.acronym(r), "_", pqname, ".txt", sep = "")
      getFunctionalAnnotationChartFile(david, fname, threshold=t, count=2L)
    }
  })
})
