setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Probe info (entrez_id)
entrez_id <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")[ , 6]

#Read expression data, assign row- and columnnames
brainNames <- list("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
brainExpr <- lapply(brainNames, function (x) {
  expr <- read.csv(paste("../ABA_human_processed/gene_expr_normalized_microarray_", x, "_2014-11-11.csv", sep = ""), header = FALSE)
  rownames(expr) <- entrez_id
  colnames(expr) <- read.csv(paste("../ABA_human_processed/sample_info_normalized_microarray_", x, "_2014-11-11.csv", sep = ""))[ , 1]
  expr
})
save(brainExpr, file = "BrainExpr.RData")

# Sample structures
ontology <- read.csv("ABA_human_processed/Ontology_edited.csv")
structures <- c("brain","frontal lobe", "mesencephalon", "striatum", "hypothalamus", "pons", "parietal lobe", "cerebellar cortex")
structureIDs <- ontology[ontology$name %in% structures, ][ , c(1:3)]

#Select PolyQ genes from data
polyQgenes <- c("HTT", "ATN1", "AR", "ATXN1", "ATXN2", "ATXN3", "CACNA1A", "ATXN7", "TBP")
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
pQgeneInfo <- probeInfo[probeInfo$gene_symbol %in% polyQgenes, ]
pQEntrezIDs <- as.character(pQgeneInfo[ , 6])
polyQgenes <- pQgeneInfo[ , 4]

#Cytoscape colors
library(RColorBrewer)
pQcolors <- brewer.pal(9, "Set1")
names(pQcolors) <- polyQgenes
save(structureIDs, pQgeneInfo, pQEntrezIDs, polyQgenes, pQcolors, file = "polyQ.RData")