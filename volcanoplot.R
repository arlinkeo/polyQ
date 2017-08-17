# Volcano plot
options(stringsAsFactors = FALSE)
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(ggplot2)
library(gridExtra)
library(grid)
# library(reshape2)
# library(plyr)
library(ggrepel)

load("resources/polyQ.RData")
load("../ABA_Rdata/BrainExpr.RData")
donorNames <- names(brainExpr)
names(donorNames) <- donorNames
load("resources/diffgenesPerBrain.RData")

# Mapping entrez IDs to gene symbols and vice versa
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector
name2entrezId <- function (x) {probeInfo$entrez_id[match(x, probeInfo$gene_symbol)]} #Input is vector

### Volcano plots of diff. expr. genes

volcano.plot <- function(tab, d, labels){
  ggplot(tab, aes(fold_change, pval_up, colour = info)) +
    geom_point(alpha = 1, size=1.5) +
    scale_colour_manual(values = c("0"="grey", "1"="red")) +
    geom_text_repel(label = labels, colour = "black", fontface = "italic",  size = 3, nudge_x = 0.2) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.title =  element_text(size = 12),
          plot.title = element_text(size = 12, face = "bold")
    ) +
    geom_vline(xintercept = 0, colour = "black") +
    geom_hline(yintercept = -log10(0.05), colour = "black") +
    labs(x = "log2 fold-change", y = "-log10 p-value") +
    ggtitle(paste("Donor ", tail(unlist(strsplit(d, split = "donor")), 1)))
}

v.plot.list <- function(donors, method = ""){
  lapply(donors, function(d){
    tab <- diffgenesPerBrain[[d]]
    # tab$fold_change <- log2(tab$MedianExprIn/tab$MedianExprOut)
    # tab$fold_change <- tab$MedianExprIn/tab$MedianExprOut
    tab$fold_change <- log2(2^(tab$MedianExprIn)/2^(tab$MedianExprOut))
    tab <- tab[, c("fold_change", method)]
    tab$info <- as.numeric(rownames(tab) %in% pQEntrezIDs)
    tab <- tab[order(tab$info),]# order of plotting
    tab$info <- as.factor(tab$info)
    tab$pval_up <- -log10(tab[, method])
    labels <- rep("", nrow(tab))
    labels[rownames(tab) %in% pQEntrezIDs] <- entrezId2Name(pQEntrezIDs)
    volcano.plot(tab, d, labels)
  })
}

png(file = "diffExprHD_up_volcanoplot3.png",900,600)

plotsUp <- v.plot.list(donorNames, method = "pval_up")
main = textGrob("Upregulated genes", gp=gpar(fontface="bold"))
grid.arrange(grobs = plotsUp, top = main, nrow = 2, ncol = 3)
dev.off()

png(file = "diffExprHD_down_volcanoplot3.png",900,600)
plotsDown <- v.plot.list(donorNames, method = "pval_down")
main = textGrob("Downregulated genes", gp=gpar(fontface="bold"))
grid.arrange(grobs = plotsDown, top = main, nrow = 2, ncol = 3) 
dev.off()
