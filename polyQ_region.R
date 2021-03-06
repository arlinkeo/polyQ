setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

# Sample structures
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")
structures <- c("frontal lobe", "mesencephalon", "striatum", "hypothalamus", "pons", "parietal lobe", "cerebellar cortex")
structureIDs <- ontology[ontology$name %in% structures, ][ , c(1:3)]
structureIDs[, 3] <- sapply(structureIDs[, 3], function(id){gsub(" ", "_", id)})
rownames(structureIDs) <- structureIDs$name

#Select PolyQ genes from data
polyQgenes <- c("HTT", "ATN1", "AR", "ATXN1", "ATXN2", "ATXN3", "CACNA1A", "ATXN7", "TBP")
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
pqrows <- which(probeInfo$gene_symbol %in% polyQgenes)
pQEntrezIDs <- as.character(probeInfo$entrez_id[pqrows])
polyQgenes <- probeInfo$gene_symbol[pqrows]

#Gephi node colors
library(RColorBrewer)
pQcolors <- brewer.pal(9, "Set1")
names(pQcolors) <- polyQgenes
save(structureIDs, pQEntrezIDs, polyQgenes, pQcolors, file = "resources/polyQ.RData")