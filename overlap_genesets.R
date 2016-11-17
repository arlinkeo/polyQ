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

genepairs <- t(combn(pQEntrezIDs, 2))
rownames(genepairs) <- apply(genepairs, 1, function(x){paste(entrezId2Name(x[1]), "-", entrezId2Name(x[2]), sep = "")})

# Function to count number of overlapping genes between two polyQ sets per region
overlap <- function(x) {
  apply(genepairs, 1, function(y){
    geneset1 <- x[[y[1]]]
    geneset2 <- x[[y[2]]]
    genes <- intersect(geneset1, geneset2)
    length(genes)
  })
}

#Load asssociations info from literature
associations <- read.csv("datatype_interactions.txt", sep = "\t", row.names = 1, comment.char = "#")
col_select <- c("SCA_total", "HD_total", "SCA_or_HD")
associations <- associations[, col_select, drop = FALSE]
associations <- associations[rownames(genepairs), ]

#Count number of overlapping genes between two polyQ sets, combine with associations info, and plot
pdf(file = "overlap_genesets.pdf", 12, 16)
for (i in 5:8) {
  f <- paste("resources/genesets_threshold0", i,"0.RData", sep = "")
  attach(f)
  table <- sapply(regionLs, overlap)
  detach(2)
  table <- cbind(table, associations)
  # Sort rows by associations
  SCA_and_HD <- which(bitwAnd(table$SCA_total, table$HD_total) == 1)
  SCA <- which(table$SCA_total == 1)
  only_SCA <- SCA[-which(SCA == SCA_and_HD)]
  HD <- which(table$HD_total == 1)
  only_HD <- HD[-which(HD == SCA_and_HD)]
  only_HD <- only_HD[order(table[only_HD, ]$SCA_total)]
  not_SCA_and_HD <- c(which(table$SCA_or_HD == 0), which(is.na(table$SCA_or_HD)))
  order <- c(only_SCA, SCA_and_HD, only_HD, not_SCA_and_HD)
  
  table <- table[order, ]
  par(mar = c(6, 10, 15, 4));
  labeledHeatmap(as.matrix((table > 0) + 0), xLabels = colnames(table), yLabels = make.italic(rownames(table)), 
                 setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), plotLegend = FALSE,
                 textMatrix = table, main = paste("Overlap between two polyQ gene sets with genes correlated >0.", i, sep = ""))
}
dev.off()