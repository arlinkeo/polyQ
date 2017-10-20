# select correlated genes based on a corr. threshold

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")

#Prepare data and functions
source("PolyQ_scripts/baseScript.R")
structureIDs <- structureIDs[!structureIDs$name %in% c("brain", "cerebellum"), ]
structureIDs <- rbind(HD_region = c(NA, "HDnetworkBD", "HD_region"), structureIDs)

#Load mean corr. data across 6 brains. and select based on threshold
structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- structureIDs$name

lapply(c(4:6), function(i){
  regionLs <- lapply(structures, function(x) {
    f <- paste("regional_coexpression/", x[3], "/meanCor_", x[2], ".RData", sep = "")
    print(f)
    attach(f)
    selection <- lapply(pQEntrezIDs, function(x){c(x, names(which(meanCor[x,] > i)))})
    names(selection) <- pQEntrezIDs
    detach(2)
    selection
  })
 save(regionLs, file = paste("resources/genesets_threshold0", i*10,"0.RData", sep = ""))
})

#Sort and plot table
colorder <- c("HD_region", "frontal_lobe", "parietal_lobe", "striatum", "hypothalamus", "mesencephalon", "cerebellar_cortex", "pons")
table.numbers <- dget("polyQ_scripts/tableNumbers.R")

plots <- lapply(c(4:6), function(i){
  f <- paste("resources/genesets_threshold0", i, "0.RData", sep = "")
  attach(f)
  rlTab <- sapply(regionLs, function(r){sapply(r, length)})-1
  detach(2)
  colnames(rlTab) <- gsub("_", " ", colnames(rlTab))
  rownames(rlTab) <- sapply(rownames(rlTab), entrezId2Name)
  plot <- table.numbers(rlTab, name = expression(atop("Co-expressed", "genes")), title = paste("r = 0.", i))
  plot
})

pdf(file = "regionLs_threshold.pdf", 10, 4)
par(mar = c(2, 10, 15, 4))
plots
dev.off()