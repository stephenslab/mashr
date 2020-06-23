## @title A wrapper function to `stop` call
labelled_stop = function(x, msg) 
  stop(paste(gsub("\\s+", " ", paste0(deparse(x))), msg), call.=F)

## @title Basic sanity check for covariance matrices
## @param X input matrix
check_covmat_basics = function(x) {
  label = substitute(x)
  if (!is.matrix(x))
    labelled_stop(label, "is not a matrix")
  if (!is.numeric(x))
    labelled_stop(label, "is not a numeric matrix")
  if (any(is.na(x)))
    labelled_stop(label, "cannot contain NA values")
  if (any(is.infinite(x)))
    labelled_stop(label, "cannot contain Inf values")
  if (any(is.nan(x)))
    labelled_stop(label, "cannot contain NaN values")
  if (nrow(x) != ncol(x))
    labelled_stop(label, "is not a square matrix")
  if (!isSymmetric(x, check.attributes = FALSE))
    labelled_stop(label, "is not a symmetric matrix")
  return(TRUE)
}

## @title check matrix for positive definitness
## @param X input matrix
check_positive_definite = function(x) {
  check_covmat_basics(x)
  tryCatch(chol(x),
    error = function(e) labelled_stop(substitute(x),
                                      "must be positive definite"))
  return(TRUE)
}
