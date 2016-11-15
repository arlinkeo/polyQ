# Co-expression distribution of polyQ genes

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("resources/polyQ.RData")
structureIDs <- structureIDs[!structureIDs$name %in% c("cerebellar nuclei","basal forebrain","globus pallidus"), ] # remove structures from list
structureIDs[, 3] <- sapply(structureIDs[, 3], function(x){gsub(" ", "_", x)})
structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- structureIDs$name
probeInfo <- read.csv("ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

# # Load data
# coexpr_dist <- lapply(structures, function(r){
#   fName <- paste("regional_coexpression/", r[3], "/meanCor_", r[2], ".RData", sep = "")
#   attach(fName)
#   mat <- meanCor[pQEntrezIDs, ]
#   detach(2)
#   mat
# })
# save(coexpr_dist, file = "resources/coexpr_dist.RData")
load("resources/coexpr_dist.RData")

# Plot co-expression distribution
pdf(file = "coexpr_dist.pdf", 16, 12)
#par(mfcol=c(9,7), oma = c(5,5,6,0), mai = c(0.2,0.2,0.5,0.2))
par(oma = c(8,8,12,2), mai = c(0.2,0.2,0.5,0.2))
layout(matrix(c(1:72), 9, 8), widths = c(2, rep(3, 7)), heights = rep(1, 9))
lapply(polyQgenes, function(pq){
  plot(0, type = "n", axes=F, xlab="", ylab="")
  mtext(make.italic(pq), 3, line = 0)
})
quantiles <- sapply(names(coexpr_dist), function(r){
  m <- coexpr_dist[[r]]
  res <- sapply(pQEntrezIDs, function(x){
    vector <- m[x, ]
    par(mai = c(0.15, 0.15, 0.15, 0.15))
    hist(vector, breaks = seq(-1,1,by=0.05), xlab = NULL, ylab = NULL, xlim = c(-1, 1), ylim = c(0, 4000), main = NULL, cex.axis = 0.8)
    q95 <- quantile(vector, probs = 0.95)
    abline(v = q95, col = "red")
    q95
  })
  mtext(gsub("_", " ", r), 3, line = 65, cex = 1.2)
  res
})
title("PolyQ co-expression distribution", outer = TRUE, cex.main = 3)
title(xlab = "Co-expression", ylab = "Frequency", outer = TRUE, cex.lab = 2.5)
dev.off()

#Print quantile info
quantiles
mean(quantiles)
median(quantiles)
sd(quantiles)