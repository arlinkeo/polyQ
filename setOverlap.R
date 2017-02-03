# Overlap between sets of genes/functional terms
setOverlap <- function(x) {
  setNames <- names(x)
  n <- length(setNames)
  mat <- matrix(0, n, n, dimnames = list(setNames, setNames))
  for (pq1 in setNames){
    for (pq2 in setNames){
      set1 <- x[[pq1]]
      set2 <- x[[pq2]]
      overlap <- intersect(set1, set2)
      mat[pq1, pq2] <- length(overlap)
    }
  }
  mat
}
