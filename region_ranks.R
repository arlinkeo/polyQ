setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
selectTop500 <- function(x) {names(head(-sort(-x), 500))} #Input is row/column of a corr. matrix

#####Load Data#####

#Load top500's for each region
regionLs <- split(structureIDs, seq(nrow(structureIDs)))
names(regionLs) <- gsub(" ", "_", structureIDs$name)
regionLs <- lapply(regionLs, function(x) {
  f <- paste("region_specific/", gsub(" ", "_", x[3]), "/polyQgenes_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  top500 <- apply(meanCor[pQEntrezIDs, ], 1, selectTop500)
  detach(2)
  top500
})

# add top 500 of whole brain to list
attach("whole_brain/meanCor.Rdata")
top500_wb <- apply(meanCor[pQEntrezIDs, ], 1, selectTop500)
detach(2)
regionLs <- c(list(top500_wb), regionLs)
names(regionLs)[1] <- "whole_brain"

######Ranking per polyQ gene#####

#Rankings of correlated genes in each region, per polyQ gene
names(pQEntrezIDs) <- pQEntrezIDs
rankings <- lapply(pQEntrezIDs, function (pqGene){
  print(pqGene)
  regTop500 <- sapply(regionLs, function(x) {
    x[ , pqGene]
  })
  top25s <- as.vector(regTop500[1:25, ])
  tab <- lapply(top25s, function(gene){
    newRow <- apply(regTop500, 2, function(col){which(col == gene)})
    newRow <- unlist(lapply(newRow, function(x){ if (length(x) == 0) NaN else x}))
  })
  tab <- t(as.data.frame(tab))
  rownames(tab) <- sapply(top25s, entrezId2Name)
  tab
})
save(rankings, file = "rankings.RData")

#Print and remove duplicates
duplicates <- lapply(rankings, function(x){rownames(x)[which(duplicated(rownames(x)))]})
names(duplicates) <- sapply(names(duplicates), entrezId2Name)
rankings <- lapply(rankings, function(x){x[unique(rownames(x)), ]})

#Visualize table
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}
colors <- colorRampPalette(c("red", "white"))(200)
plot.table <- function(x, main = ""){
  par(mar = c(2, 4, 12, 2));
  labeledHeatmap((x), xLabels = gsub("_", " ", colnames(x)), yLabels = make.italic(rownames(x)), colors = colors, #naColor = "white",
              xLabelsPosition = "top", setStdMargins = FALSE, xLabelsAdj = 0, zlim = c(0,500), cex.lab.y = 0.4, cex.lab.x = 0.8, 
              textMatrix = replace(x, which(is.na(x)), ""), cex.text = 0.4, main = main)
}

pdf(file = "polyQ_regionranks.pdf", 6, 16)
sapply(names(rankings), function(x){
  gene <- entrezId2Name(x)
  plot.table(rankings[[x]], main = substitute(paste(italic(gene), " co-expression")))})
dev.off()

#####Ranking per brain region#####
