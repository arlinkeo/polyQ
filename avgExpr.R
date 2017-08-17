#Average expression of polyQ gene across regional-specific samples and all donors

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)
library(RColorBrewer)
load("resources/BrainExpr.RData")
load("resources/sampleIDs.RData")
sampleIDs$brain <- NULL
sampleIDs$cerebellum <- NULL
source("PolyQ_scripts/baseScript.R")

#Average expression of a polyQ gene in a structure across donors and samples.
avgExpr <- sapply(sampleIDs, function(s){
  res <- sapply(donorNames, function(d){
    expr <- brainExpr[[d]]
    expr2 <- expr[pQEntrezIDs, as.logical(s[[d]])]
    apply(expr2, 1, mean) # Avg across region-specific samples per PQ per donor
  })
  apply(res, 1, mean) # Avg across region-specific samples and donors per PQ
})
rownames(avgExpr) <- sapply(rownames(avgExpr), entrezId2Name)
save(avgExpr, file = "resources/avgExpr.RData")
load("resources/avgExpr.RData")

pal <- colorRampPalette(rev(brewer.pal(3, "RdBu")))# pal = colorRampPalette(c('darkgreen', 'yellow', 'red'))
avgExprColor <- matrix(pal(100)[as.numeric(cut(avgExpr, breaks = 100))], length(polyQgenes), ncol(avgExpr))
rownames(avgExprColor) <- rownames(avgExpr)
colnames(avgExprColor) <- colnames(avgExpr)
save(avgExprColor, file = "resources/avgExprColor.RData")
load("resources/avgExprColor.RData")

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(4, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }	
}
colnames(avgExpr) <- gsub("_", " ", colnames(avgExpr))

pdf(file = "avgExpr.pdf", 8, 9)
par(oma = c(8,4,12,2), mai = c(0,0,0,0))
# layout(t(matrix(1:2)), widths = c(9, 1), heights = c(1,1))
# heatmap(avgExpr, col = avgExprColor, Colv = NA, Rowv = NA)

color.bar(pal(100), min(avgExpr), max(avgExpr), 
          ticks = format(c(mean(avgExpr), min(avgExpr), max(avgExpr), 
                    (mean(avgExpr)-min(avgExpr))/2+min(avgExpr), 
                    (max(avgExpr)-mean(avgExpr))/2+mean(avgExpr)), digits = 2))
dev.off()