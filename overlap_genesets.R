# Check for overlap in polyQ gene sets per region

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellar nuclei","basal forebrain","globus pallidus"), ] # remove structures from list
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

# Function to count number of overlapping genes between two polyQ sets per region
overlap <- function(x) {
  apply(genepairs, 1, function(y){
    geneset1 <- x[[y[1]]]
    geneset2 <- x[[y[2]]]
    genes <- intersect(geneset1, geneset2)
    length(genes)
  })
}

# Function to get significance of overlap using hypergeometric test.
hyper.test <- function(x) {
  apply(genepairs, 1, function(y){
    geneset1 <- x[[y[1]]]
    geneset2 <- x[[y[2]]]
    overlap <- length(intersect(geneset1, geneset2))
    ngs1 <- length(geneset1)
    ngs2 <- length(geneset2)
    totalGenes <- 19992
    if (overlap != 0){
    print(paste(cat(sapply(y, entrezId2Name), sep = "-"), 
                ": phyper(", overlap, " - 1, ", ngs1, ", 19992 - ", ngs1, ", ", ngs2, ", lower.tail = FALSE)", sep = ""))
    }
    phyper(overlap - 1, ngs1, totalGenes - ngs1, ngs2, lower.tail = FALSE)
  })
}

####### For different thresholds ##########
thresholds <- c("50" = "50", "60" = "60", "70" = "70", "80" = "80")

# Load sets of gene sets for different thresholds
geneSets <- lapply(thresholds, function(t) {
  f <- paste("resources/genesets_threshold0", t,".RData", sep = "")
  attach(f)
  gs <- regionLs
  detach(2)
  gs
})

# Plot number of overlapping genes and its significance
pdf(file = "overlap_genesets.pdf", 21, 28)
lapply(thresholds, function(t) {
  par(mar = c(6, 10, 15, 4))
  layout(matrix(c(1:2), 2, 1))
  # Count number of overlapping genes between two polyQ sets, combine with associations info
  table1 <- sapply(geneSets[[t]], overlap)
  table1 <- cbind(table1, associations)
  par(mai = c(0.5, 2, 3, 0.5))
  labeledHeatmap(as.matrix((table1 > 0) + 0), xLabels = colnames(table1), yLabels = make.italic(rownames(table1)),
                 setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), plotLegend = FALSE,
                 textMatrix = table1, main = paste("Overlap between two polyQ gene sets with genes correlated >0.", t, sep = ""))
  # Get significance of overlap gene sets with hypergeometric test
  table2 <- sapply(geneSets[[t]], hyper.test)
  table2a <- cbind(1 - table2, associations) # Manipulate plot function to get colors of significance right
  table2 <- apply(table2, c(1,2), function(x){format(x, digits = 2)})
  table2b <- cbind(table2, associations)
  #table2b <- apply(table2b, c(1,2), function(x){if (x < 0.01) {format(x, digits = 2, scientific = T)} else {x}})
  par(mai = c(0.5, 2, 3, 0.5));
  labeledHeatmap(table2a, xLabels = colnames(table2b), yLabels = make.italic(rownames(table2b)),
                 setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = blueWhiteRed(200)[100:200], plotLegend = FALSE,
                 textMatrix = table2b, main = paste("Significance of overlap between two polyQ gene sets with genes correlated >0.", t, sep = ""))
})
dev.off()

##### Only include associations based on interaction between causative and non-causative gene###
associations["CACNA1A−HTT", "HD_AAO"] <- 0
associations["ATN1−HTT", "HD_BMI"] <- 0
associations["ATN1−HTT", "HD_motor"] <- 0
associations["HTT-ATXN2", "HD_functional"] <- 0
associations["CACNA1A−HTT", "HD_depression"] <- 0
associations["HTT-ATXN7", "HD_depression"] <- 0
associations["HTT-ATXN1", "HD_depression"] <- 0
associations["CACNA1A−HTT", "HD_anxiety"] <- 0
associations["HTT-ATXN1", "HD_anxiety"] <- 0
associations["HTT-ATXN2", "HD_anxiety"] <- 0
associations["CACNA1A−HTT", "HD_irritability"] <- 0
associations["HTT-ATXN7", "HD_irritability"] <- 0
associations["HTT-ATXN2", "HD_irritability"] <- 0
associations["HTT-ATXN1", "HD_irritability"] <- 0
associations["ATN1-HTT", "HD_cognition"] <- 0
associations["HTT-ATXN3", "HD_cognition"] <- 0
associations["HTT-ATXN2", "HD_cognition"] <- 0

pdf(file = "overlap_genesets2.pdf", 21, 28)
lapply(thresholds, function(t) {
  par(mar = c(6, 10, 15, 4))
  layout(matrix(c(1:2), 2, 1))
  # Count number of overlapping genes between two polyQ sets, combine with associations info
  table1 <- sapply(geneSets[[t]], overlap)
  table1 <- cbind(table1, associations)
  par(mai = c(0.5, 2, 3, 0.5))
  labeledHeatmap(as.matrix((table1 > 0) + 0), xLabels = colnames(table1), yLabels = make.italic(rownames(table1)),
                 setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), plotLegend = FALSE,
                 textMatrix = table1, main = paste("Overlap between two polyQ gene sets with genes correlated >0.", t, sep = ""))
  # Get significance of overlap gene sets with hypergeometric test
  table2 <- sapply(geneSets[[t]], hyper.test)
  table2a <- cbind(1 - table2, associations)
  table2 <- apply(table2, c(1,2), function(x){format(x, digits = 2)})
  table2b <- cbind(table2, associations)
  #table2b <- apply(table2b, c(1,2), function(x){if (x < 0.01) {format(x, digits = 2, scientific = T)} else {x}})
  par(mai = c(0.5, 2, 3, 0.5));
  labeledHeatmap(table2a, xLabels = colnames(table2b), yLabels = make.italic(rownames(table2b)),
                 setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = blueWhiteRed(200)[100:200], plotLegend = FALSE,
                 textMatrix = table2b, main = paste("Significance of overlap between two polyQ gene sets with genes correlated >0.", t, sep = ""))
})
dev.off()

##### Add coexpression analysis of HD regions, only for gene sets with threshold > 0.5 ########
load("resources/genesets_threshold050.RData")
load("resources/genesets_threshold050_HDregion.RData")
regionLs <- c(regionLs, HD_region = list(selection)) # Add gene sets from HD region to list of brain structures

# Plot number of overlapping genes and its significance
pdf(file = "overlap_genesets3.pdf", 21, 28)
par(mar = c(6, 10, 15, 4))
layout(matrix(c(1:2), 2, 1))
# Count number of overlapping genes between two polyQ sets, combine with associations info
table1 <- sapply(regionLs, overlap)
table1 <- cbind(table1, associations)
par(mai = c(0.5, 2, 3, 0.5))
labeledHeatmap(as.matrix((table1 > 0) + 0), xLabels = colnames(table1), yLabels = make.italic(rownames(table1)),
               setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), plotLegend = FALSE,
               textMatrix = table1, main = "Overlap between two polyQ gene sets with genes correlated >0.50")
# Get significance of overlap gene sets with hypergeometric test
table2 <- sapply(regionLs, hyper.test)
table2a <- cbind(1 - table2, associations)
table2 <- apply(table2, c(1,2), function(x){format(x, digits = 2)})
table2b <- cbind(table2, associations)
par(mai = c(0.5, 2, 3, 0.5));
labeledHeatmap(table2a, xLabels = colnames(table2b), yLabels = make.italic(rownames(table2b)),
               setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = blueWhiteRed(200)[100:200], plotLegend = FALSE,
               textMatrix = table2b, main = "Significance of overlap between two polyQ gene sets with genes correlated >0.50")
dev.off()
