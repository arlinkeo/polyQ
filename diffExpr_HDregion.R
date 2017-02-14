#Differential expression in HD-associated region (Coppen2016)
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)
library(reshape2)

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
pValues <- apply(expr, 1, function(g){
  exprIn <- g[selection]
  exprOut <- g[!selection]
  res <- wilcox.test(exprIn, exprOut)
  res$p.value
})

# corrected p-values
cpv <- p.adjust(pValues, method = "bonferroni", length(pValues))
genes <- names(cpv)[which(cpv <0.05)]
pq_cpv <- cpv[genes[which(genes %in% pQEntrezIDs)]]
names(pq_cpv) <- sapply(names(pq_cpv), entrezId2Name)
pq_plotIdx <- which(polyQgenes %in% names(pq_cpv))
pq_cpv <- sapply(pq_cpv, function(x){format(x, digits =2, scientific = T)})

#Boxplots of pQ genes
data <- as.data.frame(t(expr[pQEntrezIDs, ]))
colnames(data) <- sapply(pQEntrezIDs, entrezId2Name)
data <- cbind(inHDregion, data)
v <- melt(data, id.vars="inHDregion")
pdf(file = "diffExpr_HDregion.pdf", 9)
a <- boxplot(value~inHDregion+variable, v, col = c("turquoise", "orange"), ylab = "Expression", xlab = "Gene", xaxt = 'n',
        at = c(1:(length(pQEntrezIDs)*3))[-seq(3, length(pQEntrezIDs)*3, 3)], ylim= c(0,12))
axis(1, at = seq(1.5,27,3), labels  = make.italic(polyQgenes), cex.axis = 0.8)
legend("topright", c("Outside","Inside"), fill = c("turquoise", "orange"))
text(pq_plotIdx*3-1.5, sapply(pq_plotIdx, function(x){idx <- c(x*2-1, x*2); max(a$out[which(a$group %in% idx)])+1}), 
     labels = paste("p = ", pq_cpv, sep = ""), cex = 0.5)
dev.off()