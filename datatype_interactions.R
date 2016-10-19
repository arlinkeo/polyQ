# Rank-sum test of co-expression interaction between polyQ genes and associations from literature.
# Associations from literature: Montcel (2014), Chen (2016), Raposo (2015), Stuitje (2015)
# Co-expression values between polyQ modules with top 25 genes

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

load("resources/polyQ.RData")
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

# Load module means (averaged co-expression between two modules with 25 genes each) of regions
regionLs <- split(structureIDs, seq(nrow(structureIDs)))
names(regionLs) <- gsub(" ", "_", structureIDs$name)
mm_list <- lapply(regionLs, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/moduleMeans_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  mm <- abs(moduleMeans)
  detach(2)
  mm
})
attach("regional_coexpression/whole_brain/moduleMeans.Rdata")
mm_wb <- abs(moduleMeans)
detach(2)
mm_list <- c(list(mm_wb), mm_list)
names(mm_list)[1] <- "whole_brain"
rm(mm_wb, regionLs)

#Load interaction info from literature
genotype_pairs <- read.csv("Genotype-based_associations.txt", sep = "\t", row.names = 1, comment.char = "#")
genepairs <- rownames(genotype_pairs)

#Convert co-expression matrices to vectors
mm <- sapply(mm_list, function(x){
  coexpr <- sapply(genepairs, function(y){
    genes <- unlist(strsplit(y, "-"))
    gene1 <-genes[1]
    gene2 <-genes[2]
    x[gene1, gene2]
  })
})
mm <- as.data.frame(mm)

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
plot.associations <- function(x, main = ""){ # x is a vector with binary associations from literature
  pvalues = ranksum(mm, x)
  print(pvalues)
  par(mar = c(65, 12, 13.5, 12.3));
  labeledHeatmap(pvalues, xLabels = colnames(pvalues), yLabels = "rank-sum p-values", xLabelsPosition = "top", setStdMargins = FALSE,
                 xLabelsAdj = 0, colors = "white",  textMatrix = round(pvalues, digits = 2), plotLegend = FALSE, 
                 main = main)
  par(mar = c(6, 12, 16, 4), new = TRUE);
  table <- cbind(mm, x)
  table <- table[order(-x), ]
  labeledHeatmap(table, xLabels = NULL, yLabels = make.italic(rownames(table)), setStdMargins = FALSE,
                 xLabelsAdj = 0, zlim = c(0,1), colors = blueWhiteRed(200)[100:200], textMatrix = round(table, digits = 2))
  par(mar = c(20, 10, 20, 8));
  cols <- c("whole_brain", gsub(" ", "_", structureIDs[, 3]))
  temp <- reshape(table, direction = "long", varying=cols, sep = "", v.names = "coexpression", timevar = "region", times = cols)
  boxplot(temp$coexpression~temp$x*temp$region, las = 2, col = c("yellow", "darkseagreen"),
          main = main,
          at = c(1:27)[-seq(3, 27, 3)], ylab = "co-expression", cex.axis = 1.5, cex.lab = 1.5)
}

#Plot table and boxplot of correlations for associated and non-associated interaction in literature
pdf(file = "datatype_interactions.pdf", 12, 16)
# par(mar = c(20, 10, 20, 8));
# boxplot(mm, las = 2, ylab = "co-expression", main = "Co-expression distribution of polyQ modules in different brain areas",
#         cex.axis = 1.5, cex.lab = 1.5)

#Associations bassed on age-at-onset
plot.associations(genotype_pairs[, 1], main = "Effect on age-at-onset, based on interaction model")
combined1 <- bitwOr(genotype_pairs[, 1], genotype_pairs[, 2])
plot.associations(combined1, main = "Effect on age-at-onset, based on interaction and single gene model")

#Associations bassed on age-at-onset and other phenotypes
combined2 <- bitwOr(genotype_pairs[, 1], genotype_pairs[, 3])
plot.associations(combined2, main = "Effect on age-at-onset and other phenotypes, based on interaction model")
combined3 <- bitwOr(genotype_pairs[, 3], genotype_pairs[, 4])
combined <- bitwOr(combined1, combined3)
plot.associations(combined, main = "Effect on age-at-onset and other phenotypes, based on interaction and single gene model")

dev.off()