# Single co-expression between polyQ genes in all regions incl. HD associated region (Coppen2016)
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs[ , 3] <- gsub(" ", "_", structureIDs[ , 3])
structureIDs <- rbind(structureIDs, c(NA, "HDnetworkBD", "HD_region"))

structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- structureIDs$name
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
# name2entrezId <- function (x) { row <- which(probeInfo$gene_symbol == x); probeInfo[row, 6]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

### Load single correlations between polyQ genes ###
sc_list <- lapply(structures, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/meanCor_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  sc <- meanCor[pQEntrezIDs, pQEntrezIDs]
  colnames(sc) <- sapply(colnames(sc), entrezId2Name)
  rownames(sc) <- sapply(rownames(sc), entrezId2Name)
  detach(2)
  sc
})
save(sc_list, file = "resources/sc_list.RData")

sort(sapply(sc_list, function(x){x[7,3]}))