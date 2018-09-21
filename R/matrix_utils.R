## @title A wrapper function to `stop` call
check_failed = function(x, msg) stop(paste(deparse(substitute(x)), msg), call.=F)

## @title Basic sanity check for covariance matrices
## @param X input matrix
check_covmat_basics = function(x) {
    if (!is.matrix(x))
        check_failed(x, "is not a matrix")
    if (!is.numeric(x))
        check_failed(x, "is not a numeric matrix")
    if (any(is.na(x)))
        check_failed(x, "cannot contain NA values")
    if (any(is.infinite(x)))
        check_failed(x, "cannot contain Inf values")
    if (any(is.nan(x)))
        check_failed(x, "cannot contain NaN values")
    if (nrow(x) != ncol(x))
        check_failed(x, "is not a square matrix")
    if (!isSymmetric(x))
        check_failed(x, "is not a symmetric matrix")
  return(TRUE)
}

## @title check matrix for positive definitness
## @param X input matrix
check_positive_definite = function(x) {
  check_covmat_basics(x)
  tryCatch(chol(x), error = function(e) check_failed(x, "must be positive definite"))
  return(TRUE)
}
