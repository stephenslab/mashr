#' @title Compute conditional likelihoods for bhat vector.
#'
#' @description Computes vector of likelihoods for bhat for each of P
#'     prior covariances.
#'
#' @param bhat bhat vector (length R)
#'
#' @param V R x R covariance matrix for likelihood.
#'
#' @param Ulist list of prior covariance matrices.
#'
#' @param log If \code{TRUE}, the return value is a matrix of
#'     log-likelihoods.
#'
#' @return Vector of length P in which the pth element contains the
#'     multivariate normal likelihood p(bhat | Ulist[[p]], V).
#'
#' @importFrom mvtnorm dmvnorm
#'
#' @export
calc_lik_vector <- function(bhat,V,Ulist,log = FALSE)
  sapply(Ulist,function(p) dmvnorm(bhat,sigma = p + V,log = log))


#' @title Compute matrix of conditional likelihoods.
#'
#' @description computes matrix of condition likelihoods for each of J
#'     rows of Bhat for each of P prior covariances.
#'
#' @param data A \code{mash} data object; e.g., created by
#'     \code{\link{set_mash_data}}.
#'
#' @param Ulist List containing the prior covariance matrices.
#'
#' @param log If \code{TRUE}, the return value is a matrix of log-
#'     likelihoods.
#'
#' @param mc.cores The argument supplied to 
#'     \code{\link[parallel]{mclapply}} specifying the number of cores
#'     to use. Note that this is only has an effect for the Rcpp version.
#'
#' @param algorithm.version Indicate R or Rcpp version
#'
#' @return J x P matrix of multivariate normal likelihoods, p(bhat |
#'     Ulist[p], V).
#'
#' @useDynLib mashr
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#'
#' @export
calc_lik_matrix <- function (data, Ulist, log = FALSE, mc.cores = 1,
                             algorithm.version = c("Rcpp","R")) {
  algorithm.version <- match.arg(algorithm.version)

  if (mc.cores > 1 & algorithm.version != "Rcpp")
    stop("Argument \"mc.cores\" only works for Rcpp version.")
  
  if (algorithm.version == "R") {
    if(is_common_cov(data)){
      res <- calc_lik_matrix_common_cov(data,Ulist,log)
      if (nrow(res) == 1)
        res <- matrix(res)
      if (ncol(res) > 1)
        colnames(res) <- names(Ulist)
    } else {
      
      # Run the (considerably slower) version that is completely
      # implemented using existing R functions.
      res <- t(sapply(1:n_effects(data),
                    function(j) calc_lik_vector(data$Bhat[j,],get_cov(data,j),
                                                Ulist,log)))
    }
    if (nrow(res) == 1)
      res <- matrix(res)
    return(res)
  }
  else if (algorithm.version == "Rcpp") {

    # Run the C implementation using the Rcpp interface.
    res <- calc_lik_rcpp(t(data$Bhat),t(data$Shat),data$V,
                         simplify2array(Ulist),log,
                         is_common_cov(data))
    res <- res$data
    
    # Get column names for R > 1.
    if (ncol(res) > 1)
      colnames(res) <- names(Ulist)
    return(res)
  }
  else
    stop("Algorithm version should be either \"R\" or \"Rcpp\"")
}

#' @title Calculate matrix of relative likelihoods.
#'
#' @description Computes matrix of relative likelihoods for each of J
#'     rows of Bhat for each of P prior covariances.
#'
#' @param data A \code{mash} data object; e.g., created by
#'     \code{\link{set_mash_data}}.
#'
#' @param Ulist List containing the prior covariance matrices.
#'
#' @param algorithm.version indicates R or Rcpp
#'
#' @return The return value is a list containing the following components:
#'
#'     \item{lik_matrix}{J x P matrix containing likelihoods p(bhat[j]
#'       + Ulist[p], V), but normalized so that the largest entry in
#'       each row is 1.}
#'
#'     \item{lfactors}{Vector which will recover the original
#'       likelihoods; for example, \code{lfactors[i] +
#'       log(lik_matrix[i,])} yields the log-likelihoods corresponding
#'       to row i of the Bhat data matrix.}
#'
#' @export
calc_relative_lik_matrix <-
  function (data, Ulist, algorithm.version= c("Rcpp","R")) {

  algorithm.version <- match.arg(algorithm.version)

  # Compute the J x P matrix of conditional log-likelihoods.
  matrix_llik <- calc_lik_matrix(data,Ulist,log = TRUE, algorithm.version)

  # Avoid numerical issues (overflow or underflow) by subtracting the
  # largest entry in each row.
  lfactors    <- apply(matrix_llik,1,max)
  matrix_llik <- matrix_llik - lfactors
  return(list(lik_matrix = exp(matrix_llik),
              lfactors   = lfactors))
}
