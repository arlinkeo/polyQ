# Check for overlap in polyQ gene sets per region
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("brain", "cerebellum"), ]
structureIDs <- rbind(c(NA, "HDregion", "HD_region"), structureIDs)
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
#entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
entrezId2Name <- function (x) { row <- match(x, probeInfo$entrez_id); probeInfo[row, 4]}
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

# ##### Only include associations based on interaction between causative and non-causative gene #####
# associations["CACNA1A−HTT", "HD_AAO"] <- 0
# associations["ATN1−HTT", "HD_BMI"] <- 0
# associations["ATN1−HTT", "HD_motor"] <- 0
# associations["HTT-ATXN2", "HD_functional"] <- 0
# associations["CACNA1A−HTT", "HD_depression"] <- 0
# associations["HTT-ATXN7", "HD_depression"] <- 0
# associations["HTT-ATXN1", "HD_depression"] <- 0
# associations["CACNA1A−HTT", "HD_anxiety"] <- 0
# associations["HTT-ATXN1", "HD_anxiety"] <- 0
# associations["HTT-ATXN2", "HD_anxiety"] <- 0
# associations["CACNA1A−HTT", "HD_irritability"] <- 0
# associations["HTT-ATXN7", "HD_irritability"] <- 0
# associations["HTT-ATXN2", "HD_irritability"] <- 0
# associations["HTT-ATXN1", "HD_irritability"] <- 0
# associations["ATN1-HTT", "HD_cognition"] <- 0
# associations["HTT-ATXN3", "HD_cognition"] <- 0
# associations["HTT-ATXN2", "HD_cognition"] <- 0

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

# Print shared co-expressed genes
ube2Fam <- read.table(file = "UBE2_genefamily.txt", header = TRUE, sep = "\t")
dnaBindGenes <- read.table(file = "dna_binding_genes.txt")[,1]
ubiqGenes <- read.table(file = "ubiquitination_genes.txt")[,1]

structs <- structureIDs$name
names(structs) <- structs

hipGene <- lapply(structs, function(r){
  sapply(polyQgenes, function(g) {
    set <- regionLs[[r]][[g]]
    set[set %in% "9026"]
  })
})
hipGene <- lapply(hipGene, function(x){x[lapply(x, length)>0]})

dnaBGenes <- lapply(structs, function(r){
  sapply(polyQgenes, function(g) {
    intersect(regionLs[[r]][[g]], dnaBindGenes)
  })
})
dnaGenesTable <- sapply(dnaBGenes, function(r){sapply(r, length)})
dnaBGenes <- lapply(dnaBGenes, function(x){x[lapply(x, length)>0]})
pdf(file = "dna_repair_genes.pdf", 8, 9)
par(mar = c(2,6,12,3));
labeledHeatmap(replace(dnaGenesTable, which(dnaGenesTable == 0), NA), xLabels = colnames(dnaGenesTable), xLabelsPosition = "top", 
                 yLabels = make.italic(rownames(dnaGenesTable)), colors = blueWhiteRed(200)[100:200], 
                 main = "DNA repair genes in gene sets", setStdMargins = FALSE, xLabelsAdj = 0, textMatrix = dnaGenesTable)
dev.off()
ubGenes <- lapply(structs, function(r){
  sapply(polyQgenes, function(g) {
    intersect(regionLs[[r]][[g]], ubGenes)
  })
})
ubGenesTable <- sapply(ubGenes, function(r){sapply(r, length)})
#dnaBGenes <- lapply(dnaBGenes, function(x){x[lapply(x, length)>0]})
pdf(file = "ubiquitination_genes.pdf", 8, 9)
par(mar = c(2,6,12,3));
labeledHeatmap(replace(dnaGenesTable, which(dnaGenesTable == 0), NA), xLabels = colnames(dnaGenesTable), xLabelsPosition = "top", 
               yLabels = make.italic(rownames(dnaGenesTable)), colors = blueWhiteRed(200)[100:200], 
               main = "Ubiquitination genes in gene sets", setStdMargins = FALSE, xLabelsAdj = 0, textMatrix = dnaGenesTable)
dev.off()

ubGenes <- sapply(polyQgenes, function(g) {
  set <- sapply(regionLs$HD_region[[g]], entrezId2Name)
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

# Plot number of overlapping genes and its significance for gene sets with threshold > 0.5
pdf(file = "overlap_geneSets.pdf", 21, 28)
par(mar = c(6, 10, 15, 4))
layout(matrix(c(1:2), 2, 1))
par(mai = c(0.5, 2, 3, 0.5))

overlapTable <- sapply(geneSetOverlap, function(r){sapply(r, length)})
overlapTable <- overlapTable[rownames(associations), ]
table1 <- cbind(overlapTable, associations)
labeledHeatmap(as.matrix((table1 > 0) + 0), xLabels = colnames(table1), yLabels = make.italic(rownames(table1)),
               setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = c("white", "red"), plotLegend = FALSE,
               textMatrix = table1, main = "Overlap between two polyQ gene sets with genes correlated >0.50")

signifTable <- geneSetOverlapSignif[rownames(associations), ]
table2a <- cbind(1 - signifTable, associations)
table2 <- apply(signifTable, c(1,2), function(x){format(x, digits = 2)})
table2b <- cbind(table2, associations)
par(mai = c(0.5, 2, 3, 0.5));
labeledHeatmap(table2a, xLabels = colnames(table2b), yLabels = make.italic(rownames(table2b)),
               setStdMargins = FALSE, xLabelsPosition = "top", xLabelsAdj = 0, colors = blueWhiteRed(200)[100:200], plotLegend = FALSE,
               textMatrix = table2b, main = "Significance of overlap between two polyQ gene sets with genes correlated >0.50")
dev.off()