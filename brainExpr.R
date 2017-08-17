#Load brain expression data
source("PolyQ_scripts/baseScript.R")

#Read expression data, assign row- and columnnames
brainExpr <- lapply(donorNames, function (x) {
  expr <- read.csv(paste("../ABA_human_processed/gene_expr_normalized_microarray_", x, "_2014-11-11.csv", sep = ""), header = FALSE)
  rownames(expr) <- probeInfo$entrez_id
  colnames(expr) <- read.csv(paste("../ABA_human_processed/sample_info_normalized_microarray_", x, "_2014-11-11.csv", sep = ""))[ , 1]
  expr
})
save(brainExpr, file = "resources/BrainExpr.RData")