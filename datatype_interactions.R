# Rank-sum test of co-expression interaction between polyQ genes and assciations from literature.
# Associations from literature: Montcel (2014), Chen (2016), Raposo (2015), Stuitje (2015)
# Co-expression values between polyQ modules with top 25 genes

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

load("resources/polyQ.RData")

# Load module means (averaged co-expression between two modules with 25 genes each) of regions
regionLs <- split(structureIDs, seq(nrow(structureIDs)))
names(regionLs) <- gsub(" ", "_", structureIDs$name)
mm_list <- lapply(regionLs, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/moduleMeans_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  mm <- moduleMeans
  #mm <- abs(moduleMeans)
  detach(2)
  mm
})
attach("regional_coexpression/whole_brain/moduleMeans.Rdata")
mm_wb <- moduleMeans
detach(2)
mm_list <- c(list(mm_wb), mm_list)
names(mm_list)[1] <- "whole_brain"
rm(mm_wb, regionLs)

#Load interaction info from literature
genotype_pairs <- read.csv("Genotype-based_interactions.txt", sep = "\t", row.names = 1, comment.char = "#")
genepairs <- rownames(genotype_pairs)
genotype_pairs <- genotype_pairs[, 1]
names(genotype_pairs) <- genepairs

#Convert co-expression matrices to vectors
mm <- sapply(mm_list, function(x){
  coexpr <- sapply(genepairs, function(y){
    genes <- unlist(strsplit(y, "-"))
    gene1 <-genes[1]
    gene2 <-genes[2]
    x[gene1, gene2]
  })
})

#Rank-sum test
res <- apply(mm, 2, function(x){
  wilcox.test(x, genotype_pairs)
})
sapply(res, function(x){x$p.value})