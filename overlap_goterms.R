# Check for overlap in functional terms of gene sets per brain region
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
# library(WGCNA)
library(reshape)
# library(ggplot2)
options(stringsAsFactors = FALSE)

#Prepare data and functions
source("PolyQ_scripts/baseScript.R")
structureIDs <- structureIDs[!structureIDs$name %in% c("brain", "cerebellum"), ]
structureIDs <- rbind(c(NA, "HDregion", "HD_region"), structureIDs)
rownames(structureIDs) <- structureIDs$name
structures <- split(structureIDs, seq(nrow(structureIDs)))
names(structures) <- structureIDs$name
setOverlap <- dget("polyQ_scripts/setOverlap.R")
setOverlapSignif <- dget("polyQ_scripts/setOverlapSignif.R")

#Function to read Rdavid output
read.RdavidOutput <- function(fileName){
  if (file.exists(fileName)){
    terms <- read.csv(fileName, header = TRUE, sep = "\t", colClasses = "character")
    if (nrow(terms) == 0){
      print("...Removed")
      file.remove(fileName)
      NULL
    } else {
      terms
    }
  } else {
    NULL
  }
}

# Load GO terms
names(pQEntrezIDs) <- pQEntrezIDs
ll <- lapply(structures, function(r){
  lapply(pQEntrezIDs, function(pq){
    pqName <- entrezId2Name(pq)
    fName <- paste("regional_coexpression/", r[3], "/goterms050_", r[2], "_", pqName, ".txt", sep = "")
    print(fName)
    read.RdavidOutput(fName)
  })
})
ll <- lapply(ll, function(s){
  names(s) <- sapply(names(s), entrezId2Name)# pQ genes to entrezID
  s
})

# Lists without multiple testing
# ll1 <- lapply(ll, function(r){lapply(r, function(pq){pq$Term})})
# Filter terms based on Benjamini p-value < 0.05
ll2 <- lapply(structures, function(id){
  n <- unlist(id[3])
  a <- unlist(id[2]) 
  r <- ll[[n]]
  dataList <- lapply(names(r), function(pq){
    dat <- r[[pq]]
    rows <- dat$Benjamini < 0.05
    res <- NULL
    if (TRUE %in% rows) {
      res <- dat[rows, c("Term", "Benjamini")]
      res <- res[order(as.numeric(res$Benjamini)), ]
      res[["Term"]] <- sapply(res$Term, function(x){unlist(strsplit(x, split = "~"))[2]})
      res[["Benjamini"]] <- sapply(res$Benjamini, function(x){format(as.numeric(x), digits =2, scientific = T)})
      write.table(res, file = paste("regional_coexpression/", n, "/signficantTerms050_", a, "_", pq, ".txt", sep = ""),
                  quote = FALSE, row.names = FALSE, sep = "\t") # file with only terms significant after multiple testing
    }
    res
  })
  names(dataList) <- names(r)
  dataList
})
table <- sapply(ll2, function(x){
  sapply(x, function(l){
    if (is.null(l)) 0 else nrow(l)
  })
})
colnames(table) <- gsub("_", " ", colnames(table))

# Plot table with number of terms in each geneset
table.numbers <- dget("polyQ_scripts/tableNumbers.R")

pdf(file = "number_of_goterms050.pdf", 10, 4)
par(mar = c(2,6,12,3));
table.numbers(table, name = expression(atop("Enriched", "functional terms")))
dev.off()

#Number of overlapping terms for each polyQ pair
ll3 <- lapply(ll2, function(r){sapply(r, function(s){s$Term})})
termSetOverlap <- lapply(ll3, setOverlap)
save(termSetOverlap, file = "resources/termSetOverlap.RData")
load("resources/termSetOverlap.RData")

#Print shared terms between pairs
lapply(structureIDs$name, function(r){
  sets <- termSetOverlap[[r]]
  pairs <- names(which(sapply(sets, length) > 10)) # select significant pairs
  sets <- sets[pairs]
  
  fileConn <- file("overlapTermSets_HDregion.txt")
  comment <- ("#PolyQ pairs with >10 shared enriched terms")
  header <- paste("PolyQ_pair", "P-value", "Overlapping_co-expressed_genes", sep = "\t")
  printList <- lapply(pairNames, function(p){
    set <- paste(sets[[p]], collapse = ", ")
    paste(p, signifPairs[[p]], set, sep = "\t")
  })
  writeLines(c(comment, header, unlist(printList)), fileConn)
  close(fileConn)
  
})

# Terms shared between HTT, ATN1, and ATXN2 in all regions
tabList <- lapply(names(ll2), function(x){
  atn1 <- ll2[[x]]$ATN1
  atxn2 <- ll2[[x]]$ATXN2
  htt <- ll2[[x]]$HTT
  terms <- Reduce(intersect, list(atn1$Term, atxn2$Term, htt$Term))
  table <- data.frame(Term = terms)
  table$atn1_p <- atn1$Benjamini[which(atn1$Term %in% terms)]
  table$atxn2_p <- atxn2$Benjamini[which(atxn2$Term %in% terms)]
  table$htt_p <- htt$Benjamini[which(htt$Term %in% terms)]
  table
})
names(tabList) <- names(ll2)
tabList <- tabList[which(sapply(tabList, nrow) != 0)] # remove empty tables


fileConn <- file("overlapTerms_HTT_ATN1_ATXN2.txt")
printList <- lapply(names(tabList), function(r){
  tab <- tabList[[r]]
  header <- paste(colnames(tab), collapse = "\t")
  rowList <- apply(tab, 1, as.list)
  c(r, header, sapply(rowList, paste, collapse = "\t"), "")
})
writeLines(unlist(printList), fileConn)
close(fileConn)
