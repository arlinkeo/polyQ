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
  #mm <- moduleMeans
  mm <- abs(moduleMeans)
  detach(2)
  mm
})
attach("regional_coexpression/whole_brain/moduleMeans.Rdata")
#mm_wb <- moduleMeans
mm_wb <- abs(moduleMeans)
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
mm <- as.data.frame(mm)

#Rank-sum test
res <- apply(mm, 2, function(x){
  x_0 <- x[which(genotype_pairs == 0)]
  x_1 <- x[which(genotype_pairs == 1)]
  wilcox.test(x_0, x_1)
})
pvalues <- t(sapply(res, function(x){x$p.value}))

#Plot table and boxplot of correlations for associated and non-associated interaction in literature
pdf(file = "datatype_interactions_abs.pdf", 12, 16)
par(mar = c(65, 12, 13.5, 12.3));
labeledHeatmap(pvalues, xLabels = colnames(pvalues), yLabels = "rank-sum p-values", xLabelsPosition = "top", setStdMargins = FALSE,
               xLabelsAdj = 0, colors = "white",  textMatrix = round(pvalues, digits = 2), plotLegend = FALSE, 
               main = "Correlations between polyQ genes across different regions")
par(mar = c(6, 12, 16, 4), new = TRUE);
table <- cbind(mm, genotype_pairs)
labeledHeatmap(table, xLabels = NULL, yLabels = make.italic(rownames(table)), setStdMargins = FALSE,
               xLabelsAdj = 0, zlim = c(0,1), colors = blueWhiteRed(200)[100:200], textMatrix = round(table, digits = 2))
par(mar = c(20, 10, 20, 8));
boxplot(mm, las = 2, ylab = "co-expression", main = "Co-expression distribution of polyQ modules in different brain areas",
        cex.axis = 1.5, cex.lab = 1.5)
par(mar = c(20, 10, 20, 8));
cols <- c("whole_brain", gsub(" ", "_", structureIDs[, 3]))
temp <- reshape(table, direction = "long", varying=cols, sep = "", v.names = "coexpression", timevar = "region", times = cols)
boxplot(temp$coexpression~temp$genotype_pairs*temp$region, las = 2, col = c("yellow", "darkseagreen"),
        main = "Co-expression distribution of associated and non-associated interaction from literature",
        at = c(1:27)[-seq(3, 27, 3)], ylab = "co-expression", cex.axis = 1.5, cex.lab = 1.5)
dev.off()