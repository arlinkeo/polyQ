# Check for overlap in functional terms of gene sets per brain region

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellar nuclei","basal forebrain","globus pallidus"), ] # remove structures from list
# structureIDs[, 3] <- sapply(structureIDs[, 3], function(x){gsub(" ", "_", x)})
structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- structureIDs$name
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
name2entrezId <- function (x) { row <- which(probeInfo$gene_symbol == x); probeInfo[row, 6]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

#Load asssociations info from literature
associations <- read.csv("datatype_interactions.txt", sep = "\t", row.names = 1, comment.char = "#")
associations <- associations[ , c(15:17, 1:14)]
# Sort rows by associations
SCA_and_HD <- which(bitwAnd(associations$SCA_total, associations$HD_total) == 1)
SCA <- which(associations$SCA_total == 1)
only_SCA <- SCA[-which(SCA == SCA_and_HD)]
HD <- which(associations$HD_total == 1)
only_HD <- HD[-which(HD == SCA_and_HD)]
only_HD <- only_HD[order(associations[only_HD, ]$SCA_total)]
not_SCA_and_HD <- c(which(associations$SCA_or_HD == 0), which(is.na(associations$SCA_or_HD)))
order <- c(only_SCA, SCA_and_HD, only_HD, not_SCA_and_HD)
associations <- associations[order, ]

# prepare gene pair matrix from rownames(associations) with entrezIds
genepairs <- t(sapply(rownames(associations), function(x){as.character(sapply(unlist(strsplit(x, "-")), name2entrezId))}))

# Function to count number of overlapping GO terms between two polyQ sets per region
overlap <- function(x) {
  apply(genepairs, 1, function(y){
    termset1 <- x[[y[1]]]
    termset2 <- x[[y[2]]]
    terms <- intersect(termset1, termset2)
    length(terms)
  })
}

#Function to read Rdavid output
read.RdavidOutput <- function(fileName){
  if (file.exists(fileName)){
    terms <- read.csv(fileName, header = TRUE, sep = "\t", colClasses = "character")
    if (nrow(terms) == 0){
      print("...Removed")
      file.remove(fileName)
      NULL
    } else {
      terms
    }
  } else {
    NULL
  }
}

# Load GO terms
names(pQEntrezIDs) <- pQEntrezIDs
ll <- lapply(structures, function(r){
  lapply(pQEntrezIDs, function(pq){
    pqName <- entrezId2Name(pq)
    fName <- paste("regional_coexpression/", gsub(" ", "_", r[3]), "/goterms050_", r[2], "_", pqName, ".txt", sep = "")
    print(fName)
    read.RdavidOutput(fName)
  })
})
ll_HDregion <- lapply(pQEntrezIDs, function(pq){
  pqName <- entrezId2Name(pq)
  fName <- paste("HD_masks_Coppen2016/goterms050_HDregions", "_", pqName, ".txt", sep = "")
  print(fName)
  read.RdavidOutput(fName)
})
ll <- c(ll, HD_region = list(ll_HDregion))
rm(ll_HDregion)
# Lists without multiple testing
ll1 <- lapply(ll, function(r){lapply(r, function(pq){pq$Term})})
# Filter terms based on Benjamini p-value < 0.05
ll2 <- lapply(ll, function(r){lapply(r, function(pq){pq[pq$Benjamini < 0.05, ]$Term})})
rm(ll)

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
  table <- cbind(table, associations)
  colnames(table) <- gsub("_", " ", colnames(table))
  labeledHeatmap(as.matrix((table > 0) + 0), xLabels = colnames(table), yLabels = make.italic(rownames(table)), 
                 setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), plotLegend = FALSE,
                 textMatrix = table, main = main)
}

pdf(file = "overlap_goterms050.pdf", 21, 28)
par(mar = c(6, 10, 15, 4));
layout(matrix(c(1:2), 2, 1))
par(mai = c(0.5, 2, 3, 0.5))
plot.overlap(ll1, main = "Overlap of GO terms between two polyQ gene sets without multiple testing")
plot.overlap(ll2, main = "Overlap of GO terms between two polyQ gene sets after multiple testing (Benjamini < 0.05)")
dev.off()