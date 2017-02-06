# Significance of overlap between sets of genes/functional terms
setOverlapSignif <- function(x, total = 19992) {
  setNames <- names(x)
  n <- length(setNames)
  mat <- matrix(0, n, n, dimnames = list(setNames, setNames))
  for (pq1 in setNames){
    for (pq2 in setNames){
      set1 <- x[[pq1]]
      set2 <- x[[pq2]]
      overlap <- length(intersect(set1, set2))
      ns1 <- length(set1)
      ns2 <- length(set2)
      if (overlap != 0){
        print(paste(cat(pq1, pq2, sep = "-"), 
                    ": phyper(", overlap, " - 1, ", ns1, ", ", total, " - ", ns1, ", ", ns2, ", lower.tail = FALSE)", sep = ""))
      }
      pVal <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
      mat[pq1, pq2] <- pVal
    }
  }
  diag(mat) <- 1
  mat
  
}
