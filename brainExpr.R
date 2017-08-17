#Load brain expression data
source("PolyQ_scripts/baseScript.R")

#Read expression data, assign row- and columnnames
# brainNames <- list("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
brainExpr <- lapply(donorNames, function (x) {
  expr <- read.csv(paste("../ABA_human_processed/gene_expr_normalized_microarray_", x, "_2014-11-11.csv", sep = ""), header = FALSE)
  rownames(expr) <- entrez_id
  colnames(expr) <- read.csv(paste("../ABA_human_processed/sample_info_normalized_microarray_", x, "_2014-11-11.csv", sep = ""))[ , 1]
  expr
})
save(brainExpr, file = "resources/BrainExpr.RData")