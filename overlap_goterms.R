# Check for overlap in GO terms of gene sets per brain region

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
#library("RDAVIDWebService")
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellar nuclei","basal forebrain","globus pallidus"), ] # remove structures from list
structureIDs[, 3] <- sapply(structureIDs[, 3], function(x){gsub(" ", "_", x)})
structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- structureIDs$name
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}
#region.acronym <- function(x) {structureIDs[structureIDs$name %in% x, ]$acronym}

#Prepare polyQ pairs
genepairs <- t(combn(pQEntrezIDs, 2))
rownames(genepairs) <- apply(genepairs, 1, function(x){paste(entrezId2Name(x[1]), "-", entrezId2Name(x[2]), sep = "")})

# Function to count number of overlapping GO terms between two polyQ sets per region
overlap <- function(x) {
  apply(genepairs, 1, function(y){
    termset1 <- x[[y[1]]]
    termset2 <- x[[y[2]]]
    terms <- intersect(termset1, termset2)
    length(terms)
  })
}

#Load asssociations info from literature
asssociations <- read.csv("datatype_interactions.txt", sep = "\t", row.names = 1, comment.char = "#")
asssociations <- asssociations[, ncol(asssociations), drop = FALSE]
asssociations <- asssociations[rownames(genepairs), ]

# Load GO terms
names(pQEntrezIDs) <- pQEntrezIDs
ll <- lapply(structures, function(r){
  lapply(pQEntrezIDs, function(pq){
    pqName <- entrezId2Name(pq)
    fName <- paste("regional_coexpression/", r[3], "/goterms070_", r[2], "_", pqName, ".txt", sep = "")
    print(fName)
    goList <- if (file.exists(fName)){read.csv(fName, header = TRUE, sep = "\t")[ , 2]} else {NULL}
  })
})
table <- sapply(ll, overlap)
table <- cbind(table, asssociations)
table <- table[order(-table[, ncol(table)]), ]
pdf(file = "overlap_goterms.pdf", 12, 16)
par(mar = c(6, 10, 15, 4));
labeledHeatmap(as.matrix((table > 0) + 0), xLabels = colnames(table), yLabels = make.italic(rownames(table)), 
               setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), plotLegend = FALSE,
               textMatrix = table, main = paste("Overlap of GO terms between two polyQ gene sets with genes correlated >0.7", sep = ""))
dev.off()