#Differential expression in HD-associated region (Coppen2016)
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)
library(reshape2)
library("RDAVIDWebService")

load("resources/polyQ.RData")
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

#Concatenate binary info of all donors
load("resources/sampleIDs.RData")
inHDregion <- do.call(c, sampleIDs[["HD_region"]]) # Binary vector for each donor
selection <- as.logical(inHDregion)
  
#Concatenate expression data of all donors
load("resources/brainExpr.RData") 
expr <- do.call(cbind, brainExpr) # 19992*3702 matrix

#Rank-sum test for each gene
pValues <- t(apply(expr, 1, function(g){
  exprIn <- g[selection]
  exprOut <- g[!selection]
  res <- wilcox.test(exprIn, exprOut)
  c(median(exprIn), median(exprOut), res$p.value)
}))
colnames(pValues) <- c("Median expr. in", "Median expr. out", "p-value")

# corrected p-values
pVal <- pValues[, 3]
corrected <- p.adjust(pVal, method = "bonferroni", length(pVal))
pValues <- as.data.frame(cbind(pValues, corrected))
save(pValues, file = "resources/diffExpr_pval.RData")
load("resources/diffExpr_pval.RData")

diffGenes <- pValues$corrected <0.01
upGenes <- Reduce("&", list(diffGenes, pValues$`Median expr. in` > pValues$`Median expr. out`))
upGenes <- pValues[upGenes, ]
upGenes <- cbind(gene_symbol = sapply(rownames(upGenes), entrezId2Name), entrez_id = rownames(upGenes), upGenes)
upGenes <- upGenes[order(upGenes$corrected), ]
#write.table(upGenes, file = "diffExpr_upregulated.txt", sep = "\t", quote = FALSE, row.names = FALSE)
downGenes <- Reduce("&", list(diffGenes, pValues$`Median expr. in` < pValues$`Median expr. out`))
downGenes <- pValues[downGenes, ]
downGenes <- cbind(gene_symbol = sapply(rownames(downGenes), entrezId2Name), entrez_id = rownames(downGenes), downGenes)
downGenes <- downGenes[order(downGenes$corrected), ]
#write.table(downGenes, file = "diffExpr_downregulated.txt", sep = "\t", quote = FALSE, row.names = FALSE)

pq_pValues <- pValues[pQEntrezIDs, ]
pq_pValues <- pq_pValues[order(pq_pValues$corrected), ]
pq_pValues$corrected <- sapply(pq_pValues$corrected, function(x){format(x, digits =2, scientific = T)})
pq_order <- rownames(pq_pValues)
pq_plotIdx <- which(pq_order %in% c(rownames(upGenes), rownames(downGenes)))
rownames(pq_pValues) <- sapply(rownames(pq_pValues), entrezId2Name)

#Boxplots of pQ genes
data <- as.data.frame(t(expr[pq_order, ]))
colnames(data) <- sapply(pq_order, entrezId2Name)
data <- cbind(inHDregion, data)
v <- melt(data, id.vars="inHDregion")
pdf(file = "diffExpr_HDregion.pdf", 9)
total <- length(pq_order)
a <- boxplot(value~inHDregion+variable, v, col = c("turquoise", "orange"), ylab = "Expression", xlab = "Gene", xaxt = 'n',
        at = c(1:(total*3))[-seq(3, total*3, 3)], ylim= c(0,12))
axis(1, at = seq(1.5, total*3,3), labels  = make.italic(colnames(data)[-1]), cex.axis = 0.8)
legend("topright", c("Outside","Inside"), fill = c("turquoise", "orange"))
positionLabels <- sapply(pq_plotIdx, function(x){idx <- c(x*2-1, x*2); max(a$out[which(a$group %in% idx)])+1})
text(pq_plotIdx*3-1.5, positionLabels, labels = paste("p = ", pq_pValues$corrected[pq_plotIdx], sep = ""), cex = 0.5)
dev.off()

#Functional enrichment of up- and down-regulated genes
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl", 
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
bg_list <- probeInfo$entrez_id
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold
setTimeOut(david, 200000)

result <- addList(david, upGenes$entrez_id, idType = "ENTREZ_GENE_ID", listName = "diffExpr_upregulated", listType = "Gene")
print(result)
setCurrentBackgroundPosition(david, 1)
getFunctionalAnnotationChartFile(david, "diffExpr_upregulated_goterms.txt", threshold=t, count=2L)

result <- addList(david, downGenes$entrez_id, idType = "ENTREZ_GENE_ID", listName = "diffExpr_downregulated", listType = "Gene")
print(result)
setCurrentBackgroundPosition(david, 1)
getFunctionalAnnotationChartFile(david, "diffExpr_downregulated_goterms.txt", threshold=t, count=2L)