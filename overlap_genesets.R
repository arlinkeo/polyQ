# Check for overlap in polyQ gene sets per region
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
library(reshape)
library(ggplot2)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("brain", "cerebellum"), ]
structureIDs <- rbind(c(NA, "HDregion", "HD_region"), structureIDs)
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
#entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
entrezId2Name <- function (x) { row <- match(x, probeInfo$entrez_id); probeInfo[row, 4]}
# name2entrezId <- function (x) { row <- which(probeInfo$gene_symbol == x); probeInfo[row, 6]} #Input is single element
# make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

# setOverlap <- dget("polyQ_scripts/setOverlap.R")
# setOverlapSignif <- dget("polyQ_scripts/setOverlapSignif.R")

load("resources/genesets_threshold050.RData")
regionLs <- lapply(regionLs, function(s){
  names(s) <- entrezId2Name(names(s))# sapply(names(s), entrezId2Name)# pQ genes to entrezID
  s <- sapply(s, entrezId2Name)
  s
})

# geneSetOverlap <- lapply(regionLs, setOverlap)
# geneSetOverlapSignif <- sapply(regionLs, function(r){setOverlapSignif(r, total = 19992)})
# save(geneSetOverlap, file = "resources/geneSetOverlap.RData")
# save(geneSetOverlapSignif, file = "resources/geneSetOverlapSignif.RData")
load("resources/geneSetOverlap.RData")
load("resources/geneSetOverlapSignif.RData")

########################################################
# Print shared co-expressed genes
ube2Fam <- read.table(file = "UBE2_genefamily.txt", header = TRUE, sep = "\t")
dnaBindGenes <- read.table(file = "dna_repair_genes.txt")[,1]
ubiqGenes <- read.table(file = "ubiquitination_genes.txt")[,1]

structs <- structureIDs$name
names(structs) <- structs

hipGene <- lapply(structs, function(r){
  sapply(polyQgenes, function(g) {
    set <- regionLs[[r]][[g]]
    set[set %in% "HIP1R"]
  })
})
hipGene <- lapply(hipGene, function(x){x[lapply(x, length)>0]})

bec1Gene <- lapply(structs, function(r){
  sapply(polyQgenes, function(g) {
    set <- regionLs[[r]][[g]]
    set[set %in% "BECN1"]
  })
})
bec1Gene <- lapply(bec1Gene, function(x){x[lapply(x, length)>0]})

dnaBGenes <- lapply(structs, function(r){
  sapply(polyQgenes, function(g) {
    intersect(regionLs[[r]][[g]], dnaBindGenes)
  })
})
intersect(dnaBGenes$HD_region$ATN1, dnaBGenes$HD_region$ATXN2)
intersect(dnaBGenes$frontal_lobe$ATN1, dnaBGenes$frontal_lobe$ATXN2)
intersect(dnaBGenes$parietal_lobe$ATN1, dnaBGenes$parietal_lobe$ATXN2)
intersect(dnaBGenes$striatum$ATN1, dnaBGenes$striatum$ATXN2)
dnaGenesTable <- sapply(dnaBGenes, function(r){sapply(r, length)})
colnames(dnaGenesTable) <- gsub("_", " ", colnames(dnaGenesTable))

ubGenes <- lapply(structs, function(r){
  sapply(polyQgenes, function(g) {
    intersect(regionLs[[r]][[g]], ubiqGenes)
  })
})
ubGenesTable <- sapply(ubGenes, function(r){sapply(r, length)})

#plot numbers
table.numbers <- dget("polyQ_scripts/tableNumbers.R")

pdf(file = "dna_repair_ubiquitination_genes.pdf", 10, 4)
par(mar = c(2,6,12,3));
table.numbers(dnaGenesTable, name = expression(atop("Co-expressed", " DNA repair genes")))
table.numbers(ubGenesTable, name = expression(atop("Co-expressed", " ubiquitination genes")))
dev.off()

ubGenes <- sapply(polyQgenes, function(g) {
  set <- regionLs$HD_region[[g]]
  res1 <- c(intersect(set, ube2Fam$Approved.Symbol), intersect(set, ube2Fam$Previous.Symbols))
  res2 <- set[grep("UB", set)]
  c(res1, res2)
})
ubGenes <- ubGenes[lapply(ubGenes, length)>0]

common_interactors <- t(sapply(polyQgenes, function(g){
  hdSet <- regionLs$HD_region[[g]]
  sapply(structureIDs$name[-8], function(r){
    common <- intersect(hdSet, regionLs[[r]][[g]])
    length(common)
  })
}))

pairs <- names(geneSetOverlap$HD_region)
common_overlap <- t(sapply(pairs, function(p){
  hdSet <- geneSetOverlap$HD_region[[p]]
  sapply(structureIDs$name[-8], function(r){
    common <- intersect(hdSet, geneSetOverlap[[r]][[p]])
    length(common)
  })
}))

sharedUbGenes <- sapply(pairs, function(p){
  set <- sapply(geneSetOverlap$HD_region[[p]], entrezId2Name)
  res1 <- c(intersect(set, ube2Fam$Approved.Symbol), intersect(set, ube2Fam$Previous.Symbols))
  res2 <- set[grep("UB", set)]
  c(res1, res2)
})
sharedUbGenes <- sharedUbGenes[lapply(sharedUbGenes, length)>0]
Reduce(intersect, sharedUbGenes)
unique(Reduce(c, sharedUbGenes))

all <-Reduce(union, regionLs$HD_region)
all <- sapply(all, entrezId2Name)
all[grep("UBE2", all)]
c(intersect(all, ube2Fam$Approved.Symbol), intersect(all, ube2Fam$Previous.Symbols))
##############################
x <- Reduce(intersect, list(regionLs$HD_region$ATN1, regionLs$HD_region$ATXN2, regionLs$HD_region$HTT))
x[grep("UB", x)]

##################
#Print overlapping gene sets
lapply(structs[1], function(r){
  
  signifPairs <- names(which(geneSetOverlapSignif[, r] < 0.05)) # select significant pairs
  signifPairs <- sort(geneSetOverlapSignif[signifPairs, r]) # sorted by significance
  signifPairs <- signif(signifPairs, digits = 2)
  pairNames <- names(signifPairs)
  sets <- geneSetOverlap[[r]][pairNames]
 
  fileConn <- file("overlapGeneSets_HDregion.txt")
  comment <- ("#PolyQ pairs with significant overlap of co-expressed genes")
  header <- paste("PolyQ_pair", "P-value", "Overlapping_co-expressed_genes", sep = "\t")
  printList <- lapply(pairNames, function(p){
    set <- paste(sets[[p]], collapse = ", ")
    paste(p, signifPairs[[p]], set, sep = "\t")
  })
  writeLines(c(comment, header, unlist(printList)), fileConn)
  close(fileConn)

})