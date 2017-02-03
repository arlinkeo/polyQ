# Overlap between sets of genes/functional terms
setOverlap <- function(x) {
  setNames <- t(combn(names(x), 2))
  rownames(setNames) <- apply(setNames, 1, function(n){paste(n[1], n[2], sep = "-")})
  overlapLs <- apply(setNames, 1, function(y){
    set1 <- x[[y[1]]]
    set2 <- x[[y[2]]]
    overlap <- intersect(set1, set2)
  })
}
