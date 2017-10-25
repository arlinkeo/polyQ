# Co-expression distribution of polyQ genes
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

#Prepare data and functions
source("PolyQ_scripts/baseScript.R")
structureIDs <- rbind(HD_region = c(NA, "HDregion", "HD_region"), structureIDs)
structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- structureIDs$name

# Load data
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
layout(matrix(c(1:90), 9, 10), widths = c(2, rep(3, 9)), heights = rep(1, 9))
lapply(polyQgenes, function(pq){
  plot(0, type = "n", axes=F, xlab="", ylab="")
  mtext(make.italic(pq), 2, line = -3.5, at = .8, cex = 1.2)
})
quantiles <- sapply(names(coexpr_dist), function(r){
  m <- coexpr_dist[[r]]
  res <- sapply(pQEntrezIDs, function(x){
    vector <- m[x, ]
    vector <- vector[names(vector) !=  x] #remove corr. with itself
    par(mai = c(0.05, 0.05, 0.05, 0.05))#c(0.15, 0.15, 0.15, 0.15))
    if (x == tail(pQEntrezIDs, 1)& r == head(names(coexpr_dist), 1)){#last row & first col
      hist(vector, breaks = seq(-1,1,by=0.05), xlim = c(-1, 1), ylim = c(0, 6000), main = NULL, cex.axis = 1)
    }else if (x == tail(pQEntrezIDs, 1)){#last row
      hist(vector, breaks = seq(-1,1,by=0.05), yaxt='n', xlim = c(-1, 1), ylim = c(0, 6000), main = NULL, cex.axis = 1)
      axis(side = 2, at = seq(0,6000,1000), labels = FALSE, tick = TRUE)
    }
    else if (r == head(names(coexpr_dist), 1)){#first col
      hist(vector, breaks = seq(-1,1,by=0.05), xaxt='n', xlim = c(-1, 1), ylim = c(0, 6000), main = NULL, cex.axis = 1)
      axis(side = 1, at = seq(-1,1,0.5), labels = FALSE, tick = TRUE)
    }
    else {
      hist(vector, breaks = seq(-1,1,by=0.05), xaxt='n', yaxt='n', xlim = c(-1, 1), ylim = c(0, 6000), main = NULL, cex.axis = 1)
      axis(side = 1, at = seq(-1,1,0.5), labels = FALSE, tick = TRUE)
      axis(side = 2, at = seq(0,6000,1000), labels = FALSE, tick = TRUE)
    }
    
    q95 <- quantile(vector, probs = 0.95)
    abline(v = q95, col = "red")
    q95
  })
  mtext(gsub("_", " ", r), 3, line = 65, cex = 1.2)
  res
})
title("PolyQ co-expression distribution", outer = TRUE, cex.main = 3)
title(xlab = expression("Pearson\'s correlation ("~italic(r)~")"), outer = TRUE, cex.lab = 2, line = 3.5)
title(ylab = "Frequency", outer = TRUE, cex.lab = 2, line = -2)
dev.off()

#Print quantile info
quantiles
mean(quantiles)
median(quantiles)
sd(quantiles)
min(quantiles)
max(quantiles)