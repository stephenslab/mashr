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
    if (sum(x == t(x)) != (nrow(x)^2))
        stop(paste(substitute(x), "is not a symmetric matrix"))
}

## @title check matrix for positive definitness
## @param X input matrix
is_positive_definite = function(x) {
  check_covmat_basics(x)
  if (det(x) <= 0)
    stop(paste(substitute(x), "must be positive definite"))
}
