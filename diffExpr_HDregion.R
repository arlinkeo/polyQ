#Differential expression in HD-associated region (Coppen2016)
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)
# library("reshape2")
# library("RDAVIDWebService")
# library("metap")

load("resources/polyQ.RData")
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
# make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}
load("resources/brainExpr.RData") 
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
load("resources/sampleIDs.RData")
sampleIDs <- sampleIDs[!names(sampleIDs) %in% c("brain", "cerebellum")]
# load("resources/polyQ.RData")
# structureIDs <- structureIDs[!structureIDs$name %in% c("brain", "cerebellum"), ]
# structureIDs <- rbind(HD_region = c(NA, "HDregion", "HD_region"), structureIDs)

inHDregion <- sampleIDs[["HD_region"]] # Binary vector for each donor
# inHDregion <- lapply(inHDregion, as.logical)

#Concatenate binary info of all donors
# inHDregion_all <- do.call(c, sampleIDs[["HD_region"]]) 
# selection_all <- as.logical(inHDregion_all)

#Concatenate expression data of all donors
# expr_all <- do.call(cbind, brainExpr) # 19992*3702 matrix

diffgenesPerBrain <- lapply(donorNames, function(b){
  expr <- brainExpr[[b]]#[1:10, ]
  inHD <- as.logical(inHDregion[[b]])
  
  # exprIn <- expr[, inHD]
  # exprOut <- expr[, !inHD]
  # genes <- rownames(expr)
  # genes <- t(sapply(genes, function(g){
  #   gIn <- unlist(exprIn[g, ])
  #   gOut <- unlist(exprOut[g, ])
  #   diffUp <- wilcox.test(gIn, gOut, alternative = "greater")
  #   diffDown <- wilcox.test(gIn, gOut, alternative = "less")
  #   c(MedianExprIn = median(gIn), MedianExprOut = median(gOut), 
  #     pval_up = diffUp$p.value, pval_down = diffDown$p.value)
  # }))

  genes <- t(apply(expr, 1, function(g){
    g <- unlist(g)
    exprIn <- g[inHD]
    exprOut <- g[!inHD]
    diffUp <- wilcox.test(exprIn, exprOut, alternative = "greater")
    diffDown <- wilcox.test(exprIn, exprOut, alternative = "less")
    c(MedianExprIn = median(exprIn), MedianExprOut = median(exprOut),
      pval_up = diffUp$p.value, pval_down = diffDown$p.value)
  }))
  genes <- as.data.frame(genes)
  # genes$pval_up_cor <- p.adjust(genes$pval_up, method = "bonferroni", nrow(genes))
  # genes$pval_down_cor <- p.adjust(genes$pval_down, method = "bonferroni", nrow(genes))
  genes
})

#Select genes p-values in all brains
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
#sapply(upGenes, entrezId2Name)
allDown <- apply(diffDown, 1, function(g) {sum(g < 0.025) >= 5})
downGenes <- names(allDown)[allDown]
#sapply(downGenes, entrezId2Name)

upGenes[which(upGenes %in% pQEntrezIDs)]
downGenes[which(downGenes %in% pQEntrezIDs)]

# combine p-values
# genes_all <- rownames(expr_all)
# combP <- sapply(genes_all, function(g){
#   pvalues <- sapply(donorNames, function(b){
#     diffgenesPerBrain[[b]][g, "corrected_p"]
#   })
#   fisherP <- sumlog(pvalues)
#   unlist(fisherP[1:3])
# })
# combP <- as.data.frame(t(combP))
# combP <- combP[order(combP$p), ]

# avgDiffExpr <- apply(simplify2array(diffgenesPerBrain), 1:2, mean)
# avgDiffExpr <- as.data.frame(avgDiffExpr)
# avgDiffExpr$combinedP <- combP
# avgDiffExpr$pVal <- NULL

# save(pValues, file = "resources/diffExpr_pval.RData")
#load("resources/diffExpr_pval.RData")

#Entrez_id's of interesting genes
diffGenes <- rownames(combP)[(combP$p < 0.001)]
higherExprIn <- t(sapply(diffGenes, function(g){# T/F across donors
  exprIn <-sapply(donorNames, function(b){ 
    diffgenesPerBrain[[b]][g, "MedianExprIn"]
  })
  exprOut <-sapply(donorNames, function(b){
    diffgenesPerBrain[[b]][g, "MedianExprOut"]
  })
  exprIn > exprOut
}))
higherInAll <- apply(higherExprIn, 1, function(x){Reduce("&", as.list(x))})# Logical vector
lowerInAll <- apply(higherExprIn, 1, function(x){Reduce("&", as.list(!x))})
upGenes <- names(which(higherInAll))
upGenes <- data.frame(gene_symbol = sapply(upGenes, entrezId2Name), entrez_id = upGenes, 
                      pvalue = sapply(combP[upGenes, "p"], function(x){format(x, digits =2, scientific = T)}))
# upGenes <- upGenes[order(upGenes$pvalue), ]
# write.table(upGenes, file = "diffExpr_upregulated.txt", sep = "\t", quote = FALSE, row.names = FALSE)
downGenes <- names(which(lowerInAll))
downGenes <- data.frame(gene_symbol = sapply(downGenes, entrezId2Name), entrez_id = downGenes, 
                        pvalue = sapply(combP[downGenes, "p"], function(x){format(x, digits =2, scientific = T)}))
# downGenes <- downGenes[order(downGenes$pvalue), ]
# write.table(downGenes, file = "diffExpr_downregulated.txt", sep = "\t", quote = FALSE, row.names = FALSE)

diffPq <- intersect(diffGenes,  pQEntrezIDs)

pq_pValues <- avgDiffExpr[pQEntrezIDs, ]
pq_pValues <- pq_pValues[order(pq_pValues$corrected), ]
pq_pValues$corrected <- sapply(pq_pValues$corrected, function(x){format(x, digits =2, scientific = T)})
pq_order <- rownames(pq_pValues)
pq_plotIdx <- which(pq_order %in% c(rownames(upGenes), rownames(downGenes)))
rownames(pq_pValues) <- sapply(rownames(pq_pValues), entrezId2Name)

#Boxplots of pQ genes
pdf(file = "diffExpr_HDregion.pdf", 9)
par(oma = c(2,2,1,1))

data <- as.data.frame(t(expr_all[pq_order, ]))
colnames(data) <- sapply(pq_order, entrezId2Name)
data <- cbind(inHDregion_all, data)
v <- melt(data, id.vars="inHDregion_all")
total <- length(pq_order)
a <- boxplot(value~inHDregion_all+variable, v, col = c("turquoise", "orange"), ylab = "Expression", xlab = "Gene", xaxt = 'n',
             at = c(1:(total*3))[-seq(3, total*3, 3)], ylim= c(0,12), 
             main = "Differential expression in all brains with combined p-values")
axis(1, at = seq(1.5, total*3,3), labels  = make.italic(colnames(data)[-1]), cex.axis = 0.8)
legend("topright", c("Outside","Inside"), fill = c("turquoise", "orange"))
positionLabels <- sapply(pq_plotIdx, function(x){idx <- c(x*2-1, x*2); max(a$out[which(a$group %in% idx)])+1})
text(pq_plotIdx*3-1.5, positionLabels, labels = paste("p = ", pq_pValues$corrected[pq_plotIdx], sep = ""), cex = 0.8)

lapply(donorNames, function(b){
  expr <- brainExpr[[b]]
  dat <- as.data.frame(t(expr[pq_order, ]))
  colnames(dat) <- sapply(pq_order, entrezId2Name)
  inHD <- inHDregion[[b]]
  dat <- cbind(inHD, dat)
  v <- melt(dat, id.vars="inHD")
  total <- length(pq_order)
  a <- boxplot(value~inHD+variable, v, col = c("turquoise", "orange"), ylab = "Expression", xlab = "Gene", xaxt = 'n',
               at = c(1:(total*3))[-seq(3, total*3, 3)], ylim= c(0,12),
               main = paste("Differential expression in ", b))
  axis(1, at = seq(1.5, total*3,3), labels  = make.italic(colnames(dat)[-1]), cex.axis = 0.8)
  legend("topright", c("Outside","Inside"), fill = c("turquoise", "orange"))
})
dev.off()
# 
# #Functional enrichment of up- and down-regulated genes
# david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl", 
#                            url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
# setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
# bg_list <- probeInfo$entrez_id
# bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
# bg
# t <- 0.05 # EASE p-value threshold
# setTimeOut(david, 200000)
# 
# result <- addList(david, upGenes$entrez_id, idType = "ENTREZ_GENE_ID", listName = "diffExpr_upregulated", listType = "Gene")
# print(result)
# setCurrentBackgroundPosition(david, 1)
# getFunctionalAnnotationChartFile(david, "diffExpr_upregulated_goterms.txt", threshold=t, count=2L)
# 
# result <- addList(david, downGenes$entrez_id, idType = "ENTREZ_GENE_ID", listName = "diffExpr_downregulated", listType = "Gene")
# print(result)
# setCurrentBackgroundPosition(david, 1)
# getFunctionalAnnotationChartFile(david, "diffExpr_downregulated_goterms.txt", threshold=t, count=2L)