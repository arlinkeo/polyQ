#Differential expression in HD-associated region (Coppen2016)
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)
library("reshape2")
library("RDAVIDWebService")
# library("metap")

source("PolyQ_scripts/baseScript.R")
load("../ABA_Rdata/BrainExpr.RData")
load("resources/sampleIDs.RData")
inHDregion <- sampleIDs[["HD_region"]] # Binary vector for each donor

diffgenesPerBrain <- lapply(donorNames, function(b){
  expr <- brainExpr[[b]]#
  inHD <- as.logical(inHDregion[[b]])

  genes <- t(apply(expr, 1, function(g){
    g <- unlist(g)
    exprIn <- g[inHD]
    exprOut <- g[!inHD]
    diffUp <- wilcox.test(exprIn, exprOut, alternative = "greater")
    diffDown <- wilcox.test(exprIn, exprOut, alternative = "less")
    c(MedianExprIn = median(exprIn), MedianExprOut = median(exprOut),
      MeanExprIn = mean(exprIn), MeanExprOut = mean(exprOut),
      pval_up = diffUp$p.value, pval_down = diffDown$p.value)
  }))
  genes <- as.data.frame(genes)
  # genes$pval_up_cor <- p.adjust(genes$pval_up, method = "bonferroni", nrow(genes))
  # genes$pval_down_cor <- p.adjust(genes$pval_down, method = "bonferroni", nrow(genes))
  genes
})
save(diffgenesPerBrain, file = "resources/diffgenesPerBrain.RData")
load("resources/diffgenesPerBrain.RData")

########################################################
#Boxplots of pQ genes
total <- length(pQEntrezIDs)

pdf(file = "diffExpr_HDregion.pdf", 12, 16)
layout(matrix(c(1:6), 6, 1), heights = rep(1, 6))
par(oma = c(12,6,12,2), mai = c(0,0.6,0,0))
lapply(donorNames, function(b){
  expr <- brainExpr[[b]][pQEntrezIDs, ]
  print(paste("max: ", max(expr), ", min: ", min(expr), sep = ""))
  rownames(expr) <- sapply(rownames(expr), entrezId2Name)
  inHD <- inHDregion[[b]]
  expr <- as.data.frame(cbind(inHD, t(expr)))
  v <- melt(expr, id.vars = "inHD")
  ylab = paste("Donor", unlist(strsplit(b, split = "donor"))[2])
  a <- boxplot(value~inHD+variable, v, col = c("turquoise", "orange"), xaxt = 'n', ylab = ylab, cex.lab = 2,
               at = c(1:(total*3))[-seq(3, total*3, 3)], ylim= c(0,11), pch = '*')
  pTab <- diffgenesPerBrain[[b]]
  maxy <- max(a$out) +0.5
  lapply(c(1:9), function(i){
    pq <- pQEntrezIDs[i]
    info <- pTab[pq,]
    p <- if (info$pval_up < 0.025) info$pval_up else if (info$pval_down < 0.05) info$pval_down else ""
    p <- format(p, digits = 2, scientific = TRUE)
    text((i*3)-1.5 ,maxy, p)
  })
  if (b == donorNames[1]) legend("topright", c("Outside","Inside"), fill = c("turquoise", "orange"), cex = 2)
})
axis(1, at = seq(1.5, total*3,3), labels  = make.italic(entrezId2Name(pQEntrezIDs)), cex.axis = 2)
title("Differential expression of polyQ genes", outer = TRUE, cex.main = 2.5)
title( ylab = "Log expression", xlab = "Gene", outer = TRUE, cex.lab = 2)
dev.off()

#########################################
#Select genes p-values in all brains (genes x brain)
diffUp <- sapply(diffgenesPerBrain, function(b){
  res <- b$pval_up
  names(res) <- rownames(b)
  res
})
diffDown <- sapply(diffgenesPerBrain, function(b){
  res <- b$pval_down
  names(res) <- rownames(b)
  res
})
# select significant
allUp <- apply(diffUp, 1, function(g) {sum(g < 0.025) >= 5})# significant in 5 out of 6 brains
upGenes <- names(allUp)[allUp]
sapply(upGenes, entrezId2Name)
allDown <- apply(diffDown, 1, function(g) {sum(g < 0.025) >= 5})
downGenes <- names(allDown)[allDown]
sapply(downGenes, entrezId2Name)

#Diff. expressed polyQ genes
upGenes[which(upGenes %in% pQEntrezIDs)]
downGenes[which(downGenes %in% pQEntrezIDs)]

values <- lapply(pQEntrezIDs, function(pq){
  sapply(diffgenesPerBrain, function(b){
    b[pq,]
  })
})
names(values) <- entrezId2Name(pQEntrezIDs)
values

#Functional enrichment of up- and down-regulated genes
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
bg_list <- probeInfo$entrez_id
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold
setTimeOut(david, 200000)

result <- addList(david, upGenes, idType = "ENTREZ_GENE_ID", listName = "diffExpr_upregulated", listType = "Gene")
print(result)
setCurrentBackgroundPosition(david, 1)
getFunctionalAnnotationChartFile(david, "diffExpr_upregulated_goterms.txt", threshold=t, count=2L)
getClusterReportFile(david, "diffExpr_upregulated_termclusters.txt", type = c("Term"))

result <- addList(david, downGenes, idType = "ENTREZ_GENE_ID", listName = "diffExpr_downregulated", listType = "Gene")
print(result)
setCurrentBackgroundPosition(david, 1)
getFunctionalAnnotationChartFile(david, "diffExpr_downregulated_goterms.txt", threshold=t, count=2L)
getClusterReportFile(david, "diffExpr_downregulated_termclusters.txt", type = c("Term"))

#Benjamini-corrected GO terms
#Function to read Rdavid output
read.RdavidOutput <- function(fileName){
  if (file.exists(fileName)){
    terms <- read.csv(fileName, header = TRUE, sep = "\t", colClasses = "character")
    if (nrow(terms) == 0){
      print("...Removed")
      file.remove(fileName)
      NULL
    } else {
      terms
    }
  } else {
    NULL
  }
}

# Read annotation chart files
upTerms <- read.RdavidOutput("diffExpr_upregulated_goterms.txt")
downTerms <- read.RdavidOutput("diffExpr_downregulated_goterms.txt")
upRows <- which(upTerms$Benjamini < 0.05)
downRows <- which(downTerms$Benjamini < 0.05)

downTerms <- downTerms[downRows, c("Category", "Term", "Count", "Benjamini")]
downTerms <- downTerms[order(downTerms$Benjamini), ] # sort by Benjamini P-value
write.table(downTerms, file = "diffExpr_downregulated_corrected_goterms.txt", sep = "\t", row.names = FALSE)