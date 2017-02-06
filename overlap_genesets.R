# Check for overlap in polyQ gene sets per region

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs[, 3] <- sapply(structureIDs[, 3], function(id){gsub(" ", "_", id)})
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
name2entrezId <- function (x) { row <- which(probeInfo$gene_symbol == x); probeInfo[row, 6]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}
setOverlap <- dget("polyQ_scripts/setOverlap.R")
setOverlapSignif <- dget("polyQ_scripts/setOverlapSignif.R")

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

# Plot number of overlapping genes and its significance for different thresholds
pdf(file = "overlap_genesets.pdf", 21, 28)
lapply(thresholds, function(t) {
  par(mar = c(6, 10, 15, 4))
  layout(matrix(c(1:2), 2, 1))
  # Count number of overlapping genes between two polyQ sets, combine with associations info
  table1 <- sapply(geneSets[[t]], function(s){sapply(setOverlap(s), length)})
  rownames(table1) <- sapply(rownames(table1), function(n){paste(sapply(unlist(strsplit(n, "-")), entrezId2Name), collapse = "-")}) 
  table1 <- table1[order, ]
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

##### Only include associations based on interaction between causative and non-causative gene #####

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
  table1 <- sapply(geneSets[[t]], function(s){sapply(setOverlap(s), length)})
  rownames(table1) <- sapply(rownames(table1), function(n){paste(sapply(unlist(strsplit(n, "-")), entrezId2Name), collapse = "-")}) 
  table1 <- table1[order, ]
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

##### Add coexpression analysis of HD regions ########

load("resources/genesets_threshold050.RData")
load("resources/genesets_threshold050_HDregion.RData")
regionLs <- c(regionLs, HD_region = list(selection)) # Add gene sets from HD region to list of brain structures
# pqOrder <- polyQgenes[c(4, 3, 7, 6, 2, 1, 5, 8, 9)]
# names(pqOrder) <- pqOrder # HTT first
regionLs <- lapply(regionLs, function(s){
  names(s) <- sapply(names(s), entrezId2Name)# pQ genes to entrezID
  # s <- sapply(pqOrder, function(pq){
  #   s[[pq]]
  # })
  s
})

matLs1 <- lapply(regionLs, setOverlap)
matLs1 <- lapply(matLs1, function(m){apply(m, c(1,2), log1p)}) # natural log(1+x) of sizes
matLs2 <- lapply(regionLs, setOverlapSignif)

as.table <- function(matLs){
  sapply(matLs, function(m){
    pairs <- t(combn(rownames(m), 2))
    rownames(pairs) <- apply(pairs, 1, function(x){paste(x[1], x[2], sep = "-")})
    apply(pairs, 1, function(x){m[x[1], x[2]]})
  })
}
table1 <- as.table(matLs1)
table2 <- as.table(matLs2)

#Export to cytoscape
export.cyt <- function(m1, m2, edgeFile = ""){
  pairs <- t(combn(rownames(m1), 2))
  colnames(pairs) <- c("fromNode", "toNode")
  weight <- apply(pairs, 1, function(x){m1[x[1], x[2]]})
  pVal <- apply(pairs, 1, function(x){m2[x[1], x[2]]})
  direction <- rep("undirected", dim(pairs)[1])
  table <- cbind(pairs, weight, pVal, direction)
  colnames(table) <- c(colnames(pairs), "log(1+size)", "p-value", "direction")
  write.table(table, file = edgeFile, sep = "\t", row.names = FALSE, quote = FALSE)
}

apply(structureIDs, 1, function(id){
  id <- unlist(id)
  structName <- id[2]
  structure <- id[3]
  fName <- paste("regional_coexpression/", structure, "/pqGeneOverlapEdges_",structName,".txt", sep = "")
  export.cyt(matLs1[[structure]], matLs2[[structure]], edgeFile = fName)
})
cyt <- export.cyt(matLs1[["HD_region"]], matLs2[["HD_region"]], edgeFile = "HD_masks_Coppen2016/pqGeneOverlapEdges_HDregion.txt")

# Plot number of overlapping genes and its significance for gene sets with threshold > 0.5


pdf(file = "overlap_genesets3.pdf", 21, 28)
par(mar = c(6, 10, 15, 4))
layout(matrix(c(1:2), 2, 1))
# Count number of overlapping genes between two polyQ sets, combine with associations info
par(mai = c(0.5, 2, 3, 0.5))
labeledHeatmap(as.matrix((table1 > 0) + 0), xLabels = colnames(table1), yLabels = make.italic(rownames(table1)),
               setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), plotLegend = FALSE,
               textMatrix = table1, main = "Overlap between two polyQ gene sets with genes correlated >0.50")
# Get significance of overlap gene sets with hypergeometric test

# table2a <- cbind(1 - table2, associations)
# table2 <- apply(table2, c(1,2), function(x){format(x, digits = 2)})
# table2b <- cbind(table2, associations)
par(mai = c(0.5, 2, 3, 0.5));
labeledHeatmap(table2a, xLabels = colnames(table2b), yLabels = make.italic(rownames(table2b)),
               setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = blueWhiteRed(200)[100:200], plotLegend = FALSE,
               textMatrix = table2b, main = "Significance of overlap between two polyQ gene sets with genes correlated >0.50")
dev.off()
