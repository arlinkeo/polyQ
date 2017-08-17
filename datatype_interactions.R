# Rank-sum test of co-expression interaction between polyQ genes and associations from literature.
# Associations from literature: Montcel (2014), Chen (2016), Raposo (2015), Stuitje (2015)
# Co-expression values between polyQ modules with top 25 genes

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(WGCNA)
options(stringsAsFactors = FALSE)

source("PolyQ_scripts/baseScript.R")

# remove brain regions wit low number of samples
structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellar nuclei","basal forebrain","globus pallidus"), ]

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
mm_list <- lapply(regionLs, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/moduleMeans_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  mm <- abs(moduleMeans)
  detach(2)
  mm
})
mm <- sapply(mm_list, function(x){mat2vec(x)})
mm <- as.data.frame(mm)

load("resources/sc_list.RData")
sc_list <- sc_list[!(names(sc_list) %in% c("CbN","BF","GP"))]
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
  table <- table[complete.cases(table), ]
  table <- table[order(-table[, ncol(table)]), ]
  labeledHeatmap(table, xLabels = NULL, yLabels = make.italic(rownames(table)), setStdMargins = FALSE,
                 xLabelsAdj = 0, zlim = c(0,1), colors = blueWhiteRed(200)[100:200], textMatrix = round(table, digits = 2))
  par(mar = c(20, 10, 20, 8));
  cols <- structureIDs[,2]
  temp <- reshape(table, direction = "long", varying=cols, sep = "", v.names = "coexpression", timevar = "region", times = cols)
  boxplot(temp$coexpression~temp$x*temp$region, las = 2, col = c("yellow", "darkseagreen"),
          main = main,
          at = c(1:(length(cols)*3))[-seq(3, length(cols)*3, 3)], ylab = "co-expression", cex.axis = 1.5, cex.lab = 1.5)
}

### Plot result (direct correlations) ###
pdf(file = "datatype_interactions3.pdf", 12, 16)

#Associations bassed on age-at-onset
plot.associations(sc, genotype_pairs[, "SCA2_AAO"], main = "Effect on age-at-onset in SCA2 patients")
plot.associations(sc, genotype_pairs[, "SCA3_AAO"], main = "Effect on age-at-onset in SCA3 patients")
plot.associations(sc, genotype_pairs[, "SCA6_AAO"], main = "Effect on age-at-onset in SCA6 patients")
plot.associations(sc, genotype_pairs[, "SCA7_AAO"], main = "Effect on age-at-onset in SCA7 patients")

plot.associations(sc, genotype_pairs[, "HD_AAO"], main = "Effect on age-at-onset in HD patients")
plot.associations(sc, genotype_pairs[, "HD_BMI"], main = "Effect on BMI in HD patients")
plot.associations(sc, genotype_pairs[, "HD_motor"], main = "Effect on motor in HD patients")
plot.associations(sc, genotype_pairs[, "HD_functional"], main = "Effect on functional score in HD patients")
plot.associations(sc, genotype_pairs[, "HD_behavioral"], main = "Effect on behavioral score in HD patients")
plot.associations(sc, genotype_pairs[, "HD_depression"], main = "Effect on depression score in HD patients")
plot.associations(sc, genotype_pairs[, "HD_anxiety"], main = "Effect on anxiety score in HD patients")
plot.associations(sc, genotype_pairs[, "HD_irritability"], main = "Effect on irritability score in HD patients")
plot.associations(sc, genotype_pairs[, "HD_cognition"], main = "Effect on cognition score in HD patients")
plot.associations(sc, genotype_pairs[, "SCA_total"], main = "Effect on SCA patients")
#plot.associations(sc, genotype_pairs[, "HD_total"], main = "Effect on HD patients") # all 1's, no 0's
plot.associations(sc, genotype_pairs[, "SCA_and_HD"], main = "Effect on SCA and HD patients")

dev.off()

### Plot result (direct correlations) ###
pdf(file = "datatype_interactions4.pdf", 12, 16)

#Associations bassed on age-at-onset
plot.associations(mm, genotype_pairs[, "SCA2_AAO"], main = "Effect on age-at-onset in SCA2 patients")
plot.associations(mm, genotype_pairs[, "SCA3_AAO"], main = "Effect on age-at-onset in SCA3 patients")
plot.associations(mm, genotype_pairs[, "SCA6_AAO"], main = "Effect on age-at-onset in SCA6 patients")
plot.associations(mm, genotype_pairs[, "SCA7_AAO"], main = "Effect on age-at-onset in SCA7 patients")

plot.associations(mm, genotype_pairs[, "HD_AAO"], main = "Effect on age-at-onset in HD patients")
plot.associations(mm, genotype_pairs[, "HD_BMI"], main = "Effect on BMI in HD patients")
plot.associations(mm, genotype_pairs[, "HD_motor"], main = "Effect on motor in HD patients")
plot.associations(mm, genotype_pairs[, "HD_functional"], main = "Effect on functional score in HD patients")
plot.associations(mm, genotype_pairs[, "HD_behavioral"], main = "Effect on behavioral score in HD patients")
plot.associations(mm, genotype_pairs[, "HD_depression"], main = "Effect on depression score in HD patients")
plot.associations(mm, genotype_pairs[, "HD_anxiety"], main = "Effect on anxiety score in HD patients")
plot.associations(mm, genotype_pairs[, "HD_irritability"], main = "Effect on irritability score in HD patients")
plot.associations(mm, genotype_pairs[, "HD_cognition"], main = "Effect on cognition score in HD patients")
plot.associations(mm, genotype_pairs[, "SCA_total"], main = "Effect on SCA patients")
#plot.associations(mm, genotype_pairs[, "HD_total"], main = "Effect on HD patients") # all 1's, no 0's
plot.associations(mm, genotype_pairs[, "SCA_and_HD"], main = "Effect on SCA and HD patients")

dev.off()