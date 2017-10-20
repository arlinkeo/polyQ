# Basic functions and variables

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

###DATA INFO###

probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")
load("resources/polyQ.RData")

donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <-donorNames

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

