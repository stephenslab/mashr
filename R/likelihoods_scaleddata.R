# NOTE: All likelihood calculations in this file are for the scaled
# data data$Bhat = S^-alpha Bhat So, for example, when alpha=1 they
# are for p(Z | g) and not for p(Bhat | g)
#
# To get to likelihoods for unscaled data need to multiply by S^alpha
# As for example is done in the functions in
# likelihoods_originaldata.R However, this scaling only matters when
# comparing fits of different models - not when fitting the model
# itself (as the scaling effects all p components the same) So we only
# do it when computing log-likelihoods for public reporting, not for
# internal use.
#

#' @title Compute conditional likelihoods for bhat vector.
#'
#' @description Computes vector of likelihoods for bhat for each of P
#'   prior covariances.
#'
#'   This is an internal (non-exported) function. This help page
#'   provides additional documentation mainly intended for developers
#'   and expert users.
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
#' @keywords internal
#' 
calc_lik_vector <- function(bhat,V,Ulist,log = FALSE)
  sapply(Ulist,function(p) dmvnorm(bhat,sigma = p + V,log = log))

#' @title calc_lik_matrix_common_cov
#' 
#' @description computes matrix of likelihoods for each of J rows of
#'   Bhat for each of P prior covariances; special case when standard
#'   errors and variances are all same across j.
#' 
#'   This is an internal (non-exported) function. This help page
#'   provides additional documentation mainly intended for developers
#'   and expert users.
#' 
#' @param data a mash data object, eg as created by \code{mash_set_data}
#' 
#' @param Ulist list of prior covariance matrices
#' 
#' @param log if true computes log-likelihood
#' 
#' @return J x P vector of multivariate normal likelihoods, p(bhat |
#' Ulist[p], V), where V is same for each bhat
#' 
#' @details Compared with \code{calc_lik_matrix} this function
#' exploits fact that observations are iid in this case, so the
#' inverse covariance matrices only need to be done once, reducing
#' computation to R^3 + JR^2 instead of JR^3
#' 
#' @importFrom plyr laply
#'
#' @keywords internal
#' 
calc_lik_matrix_common_cov = function(data, Ulist, log = FALSE){
  V   <- get_cov(data,1) # all covariances are same
  res <- laply(Ulist,function(U) {
    R <- tryCatch(chol(U + V),error = function(e) FALSE)
    if (is.logical(R))
      return(rep(-Inf,nrow(data$Bhat)))
    else
      return(dmvnorm(x = data$Bhat,sigma = V + U,log = log))
   })
  dimnames(res) = NULL # just to make result identical to the non-common-cov version
  return(t(res))
}

#' @title Compute matrix of conditional likelihoods.
#'
#' @description computes matrix of condition likelihoods for each of J
#'   rows of Bhat for each of P prior covariances.
#'
#'   This is an internal (non-exported) function. This help page
#'   provides additional documentation mainly intended for developers
#'   and expert users.
#' 
#' @param data A \code{mash} data object; e.g., created by
#'     \code{\link{mash_set_data}}.
#'
#' @param Ulist List containing the prior covariance matrices.
#'
#' @param log If \code{TRUE}, the return value is a matrix of log-
#'     likelihoods.
#'
#' @param mc.cores The argument supplied to
#'     \code{openmp} specifying the number of cores
#'     to use. Note that this is only has an effect for the Rcpp version.
#'
#' @param algorithm.version Indicate R or Rcpp version
#'
#' @return J x P matrix of multivariate normal likelihoods, p(bhat |
#'     Ulist[p], V).
#'
#' @useDynLib mashr
#'
#' @importFrom Rcpp evalCpp
#'
#' @keywords internal
#' 
calc_lik_matrix <- function (data, Ulist, log = FALSE, mc.cores = 1,
                             algorithm.version = c("Rcpp","R")) {

  algorithm.version <- match.arg(algorithm.version)

  if (mc.cores > 1 & algorithm.version != "Rcpp")
    stop("Argument \"mc.cores\" only works for Rcpp version.")

  if (algorithm.version == "R") {

    # check if the rows of Shat are same
    if(data$commonV && is_common_cov_Shat(data)){
      res <- calc_lik_matrix_common_cov(data,Ulist,log)
      if (nrow(res) == 1)
        res <- matrix(res)
      if (ncol(res) > 1)
        colnames(res) <- names(Ulist)
    }else {
      # Run the (considerably slower) version that is completely
      # implemented using existing R functions.
      res <- t(sapply(1:n_effects(data),
                    function(j) calc_lik_vector(data$Bhat[j,],get_cov(data,j),
                                                Ulist,log)))
    }
    if (nrow(res) == 1)
      res <- matrix(res)

  }
  else if (algorithm.version == "Rcpp") {
    if(!data$commonV){
      stop('effect specific V has not implemented in Rcpp')
    }
    # Run the C implementation using the Rcpp interface.
    if (is.null(data$L))
        res <- calc_lik_rcpp(t(data$Bhat),t(data$Shat),data$V,
                             matrix(0,0,0), simplify2array(Ulist), 0,
                             log, is_common_cov_Shat(data), mc.cores)
    else
        res <- calc_lik_rcpp(t(data$Bhat),t(data$Shat_orig),data$V,
                             data$L, simplify2array(Ulist), 0,
                             log, is_common_cov_Shat(data), mc.cores)
    res <- res$data

    # Get column names for R > 1.
    if (ncol(res) > 1)
      colnames(res) <- names(Ulist)
  }
  else
    stop("Algorithm version should be either \"R\" or \"Rcpp\"")

  # Give a warning if any columns have -Inf likelihoods.
  rows <- which(apply(res,2,function (x) any(is.infinite(x))))
  if (length(rows) > 0)
    warning(paste("Some mixture components result in non-finite likelihoods,",
                  "either\n","due to numerical underflow/overflow,",
                  "or due to invalid covariance matrices",
                  paste(rows,collapse=", "),
                  "\n"))

  return(res)

}

#' @title Calculate matrix of relative likelihoods.
#'
#' @description Computes matrix of relative likelihoods for each of J
#'   rows of Bhat for each of P prior covariances.
#'
#'   This is an internal (non-exported) function. This help page
#'   provides additional documentation mainly intended for developers
#'   and expert users.
#' 
#' @param data A \code{mash} data object; e.g., created by
#'     \code{\link{mash_set_data}}.
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
#' @keywords internal
#' 
calc_relative_lik_matrix <-
  function (data, Ulist, algorithm.version= c("Rcpp","R")) {

  algorithm.version <- match.arg(algorithm.version)

  # Compute the J x P matrix of conditional log-likelihoods.
  matrix_llik <- calc_lik_matrix(data,Ulist,log = TRUE,
                                 algorithm.version = algorithm.version)

  # Avoid numerical issues (overflow or underflow) by subtracting the
  # largest entry in each row.
  lfactors    <- apply(matrix_llik,1,max)
  matrix_llik <- matrix_llik - lfactors
  return(list(loglik_matrix = matrix_llik,
              lfactors   = lfactors))
}
