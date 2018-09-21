## @title Basic sanity check for covariance matrices
## @param X input matrix
check_covmat_basics = function(x) {
    if (!is.matrix(x))
        stop(paste(substitute(x), "is not a matrix"))
    if (!is.numeric(x))
        stop(paste(substitute(x), "is not a numeric matrix"))
    if (any(is.na(x)))
        stop(paste(substitute(x), "cannot contain NA values"))
    if (any(is.infinite(x)))
        stop(paste(substitute(x), "cannot contain Inf values"))
    if (any(is.nan(x)))
        stop(paste(substitute(x), "cannot contain NaN values"))
    if (nrow(x) != ncol(x))
        stop(paste(substitute(x), "is not a square matrix"))
    if (!isSymmetric(x))
        stop(paste(substitute(x), "is not a symmetric matrix"))
  return(TRUE)
}

## @title check matrix for positive definitness
## @param X input matrix
check_positive_definite = function(x) {
  check_covmat_basics(x)
  tryCatch(chol(x), error = function(e) stop(paste(substitute(x), "must be positive definite")))
  return(TRUE)
}
