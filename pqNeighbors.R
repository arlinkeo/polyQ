#Search for polyQ genes in other polyQ gene sets

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
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

#Reduce gene sets to polyQ genes only
reduce2pq <- function(r){lapply(r, function(x){intersect(x, pQEntrezIDs)})}
pqInteractors <- lapply(regionLs, reduce2pq)
pqInteractors <- lapply(pqInteractors, function(r){r[which(lengths(r) != 1)]}) # remove vectors of size 1
# pqInteractors <- pqInteractors[which(lengths(pqInteractors) != 0)] # remove empty lists (regions)
pqInteractors <- lapply(pqInteractors, function(x){names <- lapply(names(x), entrezId2Name); names(x) <- names; x})
pqInteractors <- lapply(pqInteractors, function(x){lapply(x, function(y){sapply(y, entrezId2Name)})})

#Diagonal, binary matrix with 1 indicating 2 polyQ genes are neighbours (directly co-express >0.5) 
names(polyQgenes) <- polyQgenes
pqNeighbors <- lapply(pqInteractors, function(r){
  sapply(polyQgenes, function(pq1){
    interacters <- names(r)
    row <- sapply(polyQgenes, function(pq2){
      if (is.element(pq1, interacters)) {
        if (is.element(pq2, r[[pq1]][-1])) 1 else 0
      } 
      else 0
    })
    row
  })
})
save(pqNeighbors, file = "resources/pqNeighbors.RData")

#Save pdf of matrices
pdf(file = "polyQ_neighbors_matrices.pdf", 8, 9)
par(mar = c(5,6,12,3));
lapply(names(matrixList), function(r){
  mat <- matrixList[[r]]
  labeledHeatmap(mat, xLabels = make.italic(colnames(mat)), colors = c("white", "red"), plotLegend = FALSE,
                 setStdMargins = FALSE, main = paste("PolyQ neighbours in ", gsub("_", " ", r), sep =""), 
                 cex.lab = 1.3, textMatrix = mat)
})
dev.off()

####Save pdf of table with all polyQ pairs

table <- sapply(matrixList, function(mat){
  apply(genepairs, 1, function(c){
    mat[c[1], c[2]]
  })
})

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

table <- cbind(table, associations)

# prepare gene pair matrix from rownames(associations) with entrezIds
genepairs <- t(sapply(rownames(associations), function(x){as.character(unlist(strsplit(x, "-")))}))

pdf(file = "polyQ_neighbors_table.pdf", 21, 14)
par(mar = c(6, 10, 15, 4))
par(mai = c(0.5, 2, 3, 0.5))
labeledHeatmap(table, xLabels = colnames(table), yLabels = make.italic(rownames(table)),
               setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), plotLegend = FALSE,
               textMatrix = table, main = "PolyQ neighbours co-expressed >0.5")

dev.off()