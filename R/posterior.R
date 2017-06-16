#' @title posterior_cov
#' @param Vinv R x R inverse covariance matrix for the likelihood
#' @param U R x R prior covariance matrix
#' @return R x R posterior covariance matrix
#' @description If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This function returns U1.
#' @export
posterior_cov <- function(Vinv, U){
  return(U %*% solve(Vinv %*% U + diag(nrow(U))))
}

#' @title posterior_mean
#' @param bhat R vector of observations
#' @param Vinv R x R inverse covariance matrix for the likelihood
#' @param U1 R x R posterior covariance matrix, computed using posterior_cov
#' @return R vector of posterior mean
#' @description If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This function returns mu1.
#' @export
posterior_mean <- function(bhat, Vinv, U1){
  return(U1 %*% Vinv %*% bhat)
}


#' @title Compute posterior matrices.
#'
#' @description More detailed description of function goes here.
#'
#' @param data A \code{mash} data object; e.g., created by
#'     \code{\link{set_mash_data}}.
#'
#' @param Ulist List containing the prior covariance matrices.
#'
#' @param posterior_weights Vector containing the posterior
#'     probability of each mixture component in Ulist for the data
#'
#' @param algorithm.version Indicates whether to use R or Rcpp version
#'
#' @return The return value is a list containing the following
#'    components:
#'
#'    \item{PosteriorMean}{J x R matrix of posterior means.}
#'
#'    \item{PosteriorSD}{J x R matrix of posterior (marginal) standard
#'    deviations.}
#'
#'    \item{NegativeProb}{J x R matrix of posterior (marginal)
#'     probability of being negative.}
#'
#'    \item{ZeroProb}{J x R matrix of posterior (marginal) probability
#'     of being zero.}
#'
#'    \item{lfsr}{J x R matrix of local false sign rates.}
#'
#' @useDynLib mashr
#'
#' @importFrom ashr compute_lfsr
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#'
#' @export
compute_posterior_matrices <-
  function (data, Ulist, posterior_weights,
            algorithm.version = c("Rcpp","R")) {
  algorithm.version <- match.arg(algorithm.version)

  if (algorithm.version == "R") {

    compute_posterior_matrices_R_lowmem(data, Ulist, posterior_weights)

  } else if (algorithm.version == "Rcpp") {

    # Run the C implementation using the Rcpp interface.
    res  <- calc_post_rcpp(t(data$Bhat),t(data$Shat),data$V,
                           simplify2array(Ulist),t(posterior_weights))
    lfsr <- compute_lfsr(res$post_neg,res$post_zero)
    return(list(PosteriorMean = res$post_mean,
                PosteriorSD   = res$post_sd,
                lfdr          = res$post_zero,
                NegativeProb  = res$post_neg,
                lfsr          = lfsr))
  } else
    stop("Algorithm version should be either \"R\" or \"Rcpp\"")
}

#' @title compute posterior probabilities
#' @description computes posterior probabilities that each effect came from each component
#' @param pi a K vector of mixture proportions
#' @param lik_mat a JxK matrix of likelihoods
#' @return a JxK matrix of posterior probabilities, the jth row contains posteriors for jth effect
compute_posterior_weights <- function(pi, lik_mat) {
  d    <- t(pi * t(lik_mat))
  norm <- rowSums(d) # normalize probabilities to sum to 1
  return(d/norm)
}


