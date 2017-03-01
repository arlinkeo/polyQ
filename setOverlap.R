# Overlap between sets of genes/functional terms
setOverlap <- function(x) {
  setNames <- names(x)
  n <- length(setNames)
  pairs <- t(combn(setNames, 2))
  rownames(pairs) <- apply(pairs, 1, function(x){paste(x[1], "-", x[2], sep = "")})
  apply(pairs, 1, function(pq){
    set1 <- x[[pq[1]]]
    set2 <- x[[pq[2]]]
    intersect(set1, set2)
  })
}