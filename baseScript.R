# Functions

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

###DATA INFO###

probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")
# load("resources/polyQ.RData")

donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")

# Sample structures
structures <- c("brain","frontal lobe", "mesencephalon", "striatum", "hypothalamus", "pons", "parietal lobe", "cerebellar cortex", "cerebellum")
structureIDs <- ontology[ontology$name %in% structures, ][ , c(1:3)]
structureIDs[, 3] <- sapply(structureIDs[, 3], function(id){gsub(" ", "_", id)})
rownames(structureIDs) <- structureIDs$name

#Select PolyQ genes from data
polyQgenes <- c("HTT", "ATN1", "AR", "ATXN1", "ATXN2", "ATXN3", "CACNA1A", "ATXN7", "TBP")
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
pqrows <- which(probeInfo$gene_symbol %in% polyQgenes)
pQEntrezIDs <- as.character(probeInfo$entrez_id[pqrows])
polyQgenes <- probeInfo$gene_symbol[pqrows]

#PolyQ colors
library(RColorBrewer)
pQcolors <- brewer.pal(9, "Set1")
names(pQcolors) <- polyQgenes

###FUNCTIONS###

entrezId2Name <- function (x) { #Input is vector or single element
  probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]
}

name2entrezId <- function (x) { #Input is vector or single element
  probeInfo$entrez_id[match(x, probeInfo$gene_symbol)]
}


make.italic <- function(x) { #For genes
  as.expression(lapply(x, function(x) bquote(italic(.(x)))))
}

