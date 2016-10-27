# Rank-sum test of co-expression interaction between polyQ genes and associations from literature.
# Associations from literature: Montcel (2014), Chen (2016), Raposo (2015), Stuitje (2015)
# Co-expression values between polyQ modules with top 25 genes

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

load("resources/polyQ.RData")
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

#Load interaction info from literature
genotype_pairs <- read.csv("datatype_interactions.txt", sep = "\t", row.names = 1, comment.char = "#")
genepairs <- rownames(genotype_pairs)

#Function to convert co-expression matrices to vectors
mat2vec <- function(x){
  coexpr <- sapply(genepairs, function(y){
    genes <- unlist(strsplit(y, "-"))
    gene1 <-genes[1]
    gene2 <-genes[2]
    x[gene1, gene2]
  })
}

# Load module means (averaged co-expression between two modules with 25 genes each) of regions
regionLs <- split(structureIDs, seq(nrow(structureIDs)))
names(regionLs) <- gsub(" ", "_", structureIDs$acronym)
# mm_list <- lapply(regionLs, function(x) {
#   f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/moduleMeans_", x[2], ".RData", sep = "")
#   print(f)
#   attach(f)
#   mm <- abs(moduleMeans)
#   detach(2)
#   mm
# })
# mm <- sapply(mm_list, function(x){mat2vec(x)})
# mm <- as.data.frame(mm)

### Load single correlations between polyQ genes ###
sc_list <- lapply(regionLs, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/meanCor_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  sc <- abs(meanCor[pQEntrezIDs, pQEntrezIDs])
  colnames(sc) <- pQgeneInfo[pQgeneInfo$entrez_id %in% colnames(sc), "gene_symbol"]
  rownames(sc) <- pQgeneInfo[pQgeneInfo$entrez_id %in% rownames(sc), "gene_symbol"]
  detach(2)
  sc
})
rm(regionLs)
sc <- sapply(sc_list, function(x){mat2vec(x)})
sc <- as.data.frame(sc)

#Rank-sum test function
ranksum <- function(a, b) {
  res <- apply(a, 2, function(x){
    x_0 <- x[which(b == 0)]
    x_1 <- x[which(b == 1)]
    wilcox.test(x_0, x_1)
  })
  t(sapply(res, function(x){x$p.value})) # list of p-values
}

# Function to plot table and barplots
plot.associations <- function(mat, x, main = ""){ # x is a vector with binary associations from literature
  pvalues = ranksum(mat, x)
  print(pvalues)
  par(mar = c(65, 12, 13.5, 12.3));
  labeledHeatmap(pvalues, xLabels = colnames(pvalues), yLabels = "rank-sum p-values", xLabelsPosition = "top", setStdMargins = FALSE,
                 xLabelsAdj = 0, colors = "white",  textMatrix = round(pvalues, digits = 2), plotLegend = FALSE, 
                 main = main)
  par(mar = c(6, 12, 16, 4), new = TRUE);
  table <- cbind(mat, x)
  table <- table[order(-x), ]
  labeledHeatmap(table, xLabels = NULL, yLabels = make.italic(rownames(table)), setStdMargins = FALSE,
                 xLabelsAdj = 0, zlim = c(0,1), colors = blueWhiteRed(200)[100:200], textMatrix = round(table, digits = 2))
  par(mar = c(20, 10, 20, 8));
  cols <- c("Wb", gsub(" ", "_", structureIDs[, 2]))
  temp <- reshape(table, direction = "long", varying=cols, sep = "", v.names = "coexpression", timevar = "region", times = cols)
  boxplot(temp$coexpression~temp$x*temp$region, las = 2, col = c("yellow", "darkseagreen"),
          main = main,
          at = c(1:30)[-seq(3, 30, 3)], ylab = "co-expression", cex.axis = 1.5, cex.lab = 1.5)
}

### Plot results using averaged correlation between two polyQ modules ###
pdf(file = "datatype_interactions2.pdf", 12, 16)
# par(mar = c(20, 10, 20, 8));
# boxplot(mm, las = 2, ylab = "co-expression", main = "Co-expression distribution of polyQ modules in different brain areas",
#         cex.axis = 1.5, cex.lab = 1.5)

#Associations based on age-at-onset
plot.associations(mm, genotype_pairs[, 1], main = "Effect on age-at-onset, based on interaction model")
combined1 <- bitwOr(genotype_pairs[, 1], genotype_pairs[, 2])
plot.associations(mm, combined1, main = "Effect on age-at-onset, based on interaction and single gene model")

#Associations based on age-at-onset and other phenotypes
combined2 <- bitwOr(genotype_pairs[, 1], genotype_pairs[, 3])
plot.associations(mm, combined2, main = "Effect on age-at-onset and other phenotypes, based on interaction model")
combined3 <- bitwOr(genotype_pairs[, 3], genotype_pairs[, 4])
combined <- bitwOr(combined1, combined3)
plot.associations(mm, combined, main = "Effect on age-at-onset and other phenotypes, based on interaction and single gene model")

#Associations based on age-at-onset and other phenotypes in HD patients (Stuitje report)
plot.associations(mm, genotype_pairs[, 5], main = "Effect on age-at-onset and other phenotypes in HD patients, based on interaction and single gene model")

#Associations based on age-at-onset and other phenotypes in SCA patients
plot.associations(mm, genotype_pairs[, 6], main = "Effect on age-at-onset and other phenotypes in SCA patients, based on interaction and single gene model")

dev.off()
# 
# ### Plot results using single correlation between two polyQ genes ###
# pdf(file = "datatype_interactions_single.pdf", 12, 16)
# 
# #Associations bassed on age-at-onset
# plot.associations(sc, genotype_pairs[, 1], main = "Effect on age-at-onset, based on interaction model")
# combined1 <- bitwOr(genotype_pairs[, 1], genotype_pairs[, 2])
# plot.associations(sc, combined1, main = "Effect on age-at-onset, based on interaction and single gene model")
# 
# #Associations bassed on age-at-onset and other phenotypes
# combined2 <- bitwOr(genotype_pairs[, 1], genotype_pairs[, 3])
# plot.associations(sc, combined2, main = "Effect on age-at-onset and other phenotypes, based on interaction model")
# combined3 <- bitwOr(genotype_pairs[, 3], genotype_pairs[, 4])
# combined <- bitwOr(combined1, combined3)
# plot.associations(sc, combined, main = "Effect on age-at-onset and other phenotypes, based on interaction and single gene model")
# 
# dev.off()