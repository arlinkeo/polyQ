# Function to get correlation between al genes
corMatrix <- function(v, expr){ # input binary vector & expr. matrix
  samples <- as.logical(v)# row/col numbers to selec
  corMat <- cor(t(expr[ , samples]))
  diag(corMat) <- 0
  corMat
}