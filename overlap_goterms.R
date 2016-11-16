# Check for overlap in GO terms of gene sets per brain region

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
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
    fName <- paste("regional_coexpression/", r[3], "/goterms050_", r[2], "_", pqName, ".txt", sep = "")
    print(fName)
    goList <- if (file.exists(fName)){
      terms <- read.csv(fName, header = TRUE, sep = "\t", colClasses = "character")
      if (nrow(terms) == 0){
        print("...Removed")
        #file.remove(fName)
        NULL
      } else {
        terms
      }
    } else {
      NULL
      }
  })
})

# Lists without multiple testing
ll1 <- lapply(ll, function(r){lapply(r, function(pq){pq$Term})})
# Filter terms based on Benjamini p-value < 0.05
ll2 <- lapply(ll, function(r){lapply(r, function(pq){pq[pq$Benjamini < 0.05, ]$Term})})

rm(ll)# data overload

# Plot table with number of terms in each geneset
plot.numbers <- function(l, main = ""){
  table <- sapply(l, function(r){sapply(r, length)})
  rownames(table) <- sapply(rownames(table), entrezId2Name)
  Total <- apply(table, 2, sum)
  table <- rbind(table, Total)
  labeledHeatmap(replace(table, which(table == 0), NA), xLabels = gsub("_", " ", colnames(table)), xLabelsPosition = "top", 
                 yLabels = c(make.italic(rownames(table)[-10]), rownames(table)[10]), colors = blueWhiteRed(200)[100:200], 
                 main = main, setStdMargins = FALSE, xLabelsAdj = 0, textMatrix = table)
}

pdf(file = "number_of_goterms050.pdf", 8, 9)
par(mar = c(2,6,12,3));
plot.numbers(ll1, main = paste("Number of GO terms without multiple testing", sep = ))
plot.numbers(ll2, main = paste("Number of GO terms after multiple testing (Benjamini < 0.05)", sep = ))
dev.off()

# Plot table of numbers of overlapping terms
plot.overlap <- function(l, main = ""){
  table <- sapply(l, overlap)
  table <- cbind(table, asssociations)
  table <- table[order(-table[, ncol(table)]), ]
  labeledHeatmap(as.matrix((table > 0) + 0), xLabels = colnames(table), yLabels = make.italic(rownames(table)), 
                 setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), plotLegend = FALSE,
                 textMatrix = table, main = main)
}

pdf(file = "overlap_goterms050.pdf", 12, 16)
par(mar = c(6, 10, 15, 4));
plot.overlap(ll1, main = "Overlap of GO terms between two polyQ gene sets without multiple testing")
plot.overlap(ll2, main = "Overlap of GO terms between two polyQ gene sets after multiple testing (Benjamini < 0.05)")
dev.off()