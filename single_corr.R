# Single co-expression between polyQ genes in all regions incl. HD associated region (Coppen2016)

#Prepare data and functions
source("PolyQ_scripts/baseScript.R")
structureIDs <- rbind(c(NA, "HDnetworkBD", "HD_region"), structureIDs)
structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- structureIDs$name

### Load single correlations between polyQ genes ###
sc_list <- lapply(structures, function(x) {
  f <- paste("regional_coexpression/", gsub(" ", "_", x[3]), "/meanCor_", x[2], ".RData", sep = "")
  print(f)
  attach(f)
  sc <- meanCor[pQEntrezIDs, pQEntrezIDs]
  colnames(sc) <- sapply(colnames(sc), entrezId2Name)
  rownames(sc) <- sapply(rownames(sc), entrezId2Name)
  detach(2)
  sc
})
save(sc_list, file = "resources/sc_list.RData")