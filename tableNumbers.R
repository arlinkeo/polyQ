# Plot heatmap table for polyQ genes and brain structures

library(ggplot2)

tableNumbers <- function(tab, name, title = ""){
  
  tab.m <- melt(tab)
  tab.m$X2 <- factor(tab.m$X2, levels = colnames(tab)) # keep order of cols
  tab.m$X1 <- factor(tab.m$X1, levels = rev(rownames(tab))) # keep order of rows
  
  ggplot(tab.m, aes(X2, X1)) + 
    geom_tile(aes(fill = value), colour = "black") + 
    scale_fill_gradient(low = "white", high = "steelblue", name = name) + 
    geom_text(aes(label = value)) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 12), 
          axis.text.y = element_text(face = "italic")) +
    labs(x = "Brain region", y = "PolyQ gene", title = title) +
    coord_fixed(0.5)
}

