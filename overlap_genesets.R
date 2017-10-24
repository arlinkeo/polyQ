# Check for overlap in polyQ gene sets per region
source("C:/Users/dkeo/surfdrive/polyQ_coexpression/PolyQ_scripts/baseScript.R")

#Prepare data and functions
structureIDs <- rbind(c(NA, "HDregion", "HD_region"), structureIDs)

setOverlap <- dget("polyQ_scripts/setOverlap.R")
setOverlapSignif <- dget("polyQ_scripts/setOverlapSignif.R")

load("resources/genesets_threshold050.RData")
regionLs <- lapply(regionLs, function(s){
  names(s) <- entrezId2Name(names(s))# sapply(names(s), entrezId2Name)# pQ genes to entrezID
  s <- sapply(s, entrezId2Name)
  s
})

geneSetOverlap <- lapply(regionLs, setOverlap)
geneSetOverlapSignif <- sapply(regionLs, function(r){setOverlapSignif(r, total = 19992)})
save(geneSetOverlap, file = "resources/geneSetOverlap.RData")
save(geneSetOverlapSignif, file = "resources/geneSetOverlapSignif.RData")
# load("resources/geneSetOverlap.RData")
# load("resources/geneSetOverlapSignif.RData")

########################################################
# Print shared co-expressed genes
resources <- "PolyQ_scripts/resources/"
ube2Fam <- read.table(file = paste0(resources, "UBE2_genefamily.txt"), header = TRUE, sep = "\t") # HUGO nomenclature
dnaBindGenes <- read.table(file = paste0(resources, "dna_repair_genes.txt"))[,1] # MSigDB
ubiqGenes <- read.table(file = paste0(resources, "ubiquitination_genes.txt"))[,1] # MSigDB

structs <- structureIDs$name
names(structs) <- structs

# Presence of genes of which the product beclin 1 is described to interact by Ashkenazi et al. 2017
bec1Gene <- lapply(structs, function(r){
  sapply(polyQgenes, function(g) {
    set <- regionLs[[r]][[g]]
    set[set %in% "BECN1"]
  })
})
bec1Gene <- lapply(bec1Gene, function(x){x[lapply(x, length)>0]})
bec1Gene

# Number of DNA repair genes co-expressed with a polyQ gene (MSigDB)
dnaBGenes <- lapply(structs, function(r){
  sapply(polyQgenes, function(g) {
    intersect(regionLs[[r]][[g]], dnaBindGenes)
  })
})
dnaGenesTable <- sapply(dnaBGenes, function(r){sapply(r, length)})
colnames(dnaGenesTable) <- gsub("_", " ", colnames(dnaGenesTable))
dnaGenesTable

# Number of ubiquitination genes co-expressed with a polyQ gene (MSigDB)
ubGenes <- lapply(structs, function(r){
  sapply(polyQgenes, function(g) {
    intersect(regionLs[[r]][[g]], ubiqGenes)
  })
})
ubGenesTable <- sapply(ubGenes, function(r){sapply(r, length)})
colnames(ubGenesTable) <- gsub("_", " ", colnames(ubGenesTable))
ubGenesTable

#plot numbers of DNA and ubiquitination genes (manuscript, Figure S7)
table.numbers <- dget("polyQ_scripts/tableNumbers.R")

pdf(file = "dna_repair_ubiquitination_genes.pdf", 10, 4)
par(mar = c(2,6,12,3));
table.numbers(dnaGenesTable, name = expression(atop("Co-expressed", " DNA repair genes")))
table.numbers(ubGenesTable, name = expression(atop("Co-expressed", " ubiquitination genes")))
dev.off()

# Ubiquitin conjugating enzymes genes family (UBE2) 
ubGenes <- sapply(polyQgenes, function(g) {
  set <- regionLs$HD_region[[g]]
  res1 <- c(intersect(set, ube2Fam$Approved.Symbol), intersect(set, ube2Fam$Previous.Symbols)) # Also check alternative symbols
  res1
})
ubGenes[lapply(ubGenes, length)>0] # show non-empty lists

# Ubiquitin genes co-expressing with ATN1, ATXN2, and HTT in the HD-associated region
x <- Reduce(intersect, list(regionLs$HD_region$ATN1, regionLs$HD_region$ATXN2, regionLs$HD_region$HTT))
x[grep("UB", x)]

##################
#Print overlapping gene sets (manuscript, Table S3)
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