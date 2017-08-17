#Plot overlap table

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
library(reshape)
library(ggplot2)

#Load tables
load("resources/geneSetOverlap.RData")
load("resources/geneSetOverlapSignif.RData")
load("resources/termSetOverlap.RData")

#Prepare gene table
geneTab <- sapply(geneSetOverlap, function(r){sapply(r, length)})
colnames(geneTab) <- gsub("_", " ", colnames(geneTab))
geneTab.m <- melt(geneTab)
geneTab.m$X2 <- factor(geneTab.m$X2, levels = colnames(geneTab)) # keep order of cols
geneTab.m$X1 <- factor(geneTab.m$X1, levels = rev(rownames(geneTab))) # keep order of rows

# Bold labels if signficant # of genes
signif <- as.vector(geneSetOverlapSignif < 0.05)
geneTabFace <- ifelse(signif, "bold", "plain")

#Prepare term table
termTab <- sapply(termSetOverlap, function(r){sapply(r, length)})
colnames(termTab) <- gsub("_", " ", colnames(termTab))
termTab.m <- melt(termTab)
termTab.m$X2 <- factor(termTab.m$X2, levels = colnames(termTab)) # keep order of cols
termTab.m$X1 <- factor(termTab.m$X1, levels = rev(rownames(termTab))) # keep order of rows

# Bold labels if > 10 terms
atLeast10 <- termTab.m$value >= 10
termTabFace <- ifelse(atLeast10, "bold", "plain")

pdf(file = "overlap_tables.pdf", 10, 9)
par(mar = c(6, 10, 15, 4))

#Plot overlap of shared co-expressed genes (p< 0.05)
ggplot(geneTab.m, aes(X2, X1)) + 
  geom_tile(aes(fill = value), colour = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue", name = expression(atop("Shared","co-expressed genes"))) + 
  geom_text(aes(label = value, fontface = geneTabFace)) +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 12), 
        axis.text.y = element_text(face = "italic")) +
  labs(x = "Brain region", y = "PolyQ gene pair") +
  coord_fixed(0.4)

#plot overlap of function terms heatmap with >10 terms highlighted
ggplot(termTab.m, aes(X2, X1)) + 
  geom_tile(aes(fill = value), colour = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue", name = "Functional overlap") + 
  geom_text(aes(label = value, fontface = termTabFace)) +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 12), 
        axis.text.y = element_text(face = "italic")) +
  labs(x = "Brain region", y = "PolyQ gene pair") +
  coord_fixed(0.4)

dev.off()