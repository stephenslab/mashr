#' @title posterior_cov
#' @param Vinv R x R inverse covariance matrix for the likelihood
#' @param U R x R prior covariance matrix
#' @return R x R posterior covariance matrix
#' @description If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This function returns U1.
posterior_cov <- function(Vinv, U){
  return(U %*% solve(Vinv %*% U + diag(nrow(U))))
}

#' @title posterior_mean
#' @param bhat R vector of observations
#' @param Vinv R x R inverse covariance matrix for the likelihood
#' @param U1 R x R posterior covariance matrix, computed using posterior_cov
#' @return R vector of posterior mean
#' @description If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This function returns mu1.
posterior_mean <- function(bhat, Vinv, U1){
  return(U1 %*% (Vinv %*% bhat))
}

#' @title posterior_mean_matrix
#' @param Bhat J by R matrix of observations
#' @param Vinv R x R inverse covariance matrix for the likelihood
#' @param U1 R x R posterior covariance matrix, computed using posterior_cov
#' @return R vector of posterior mean
#' @description Computes posterior mean under multivariate normal model for each row of matrix Bhat.
#' Note that if bhat is N_R(b,V) and b is N_R(0,U) then b|bhat N_R(mu1,U1).
#' This function returns a matrix with jth row equal to mu1(bhat) for bhat= Bhat[j,].
posterior_mean_matrix <- function(Bhat, Vinv, U1){
  return(Bhat %*% (Vinv %*% U1))
}



#' @title Compute posterior matrices.
#'
#' @description More detailed description of function goes here.
#'
#' @param data A \code{mash} data object; e.g., created by
#'     \code{\link{mash_set_data}} or \code{\link{mash_set_data_contrast}}.
#'
#' @param Ulist List containing the prior covariance matrices.
#'
#' @param posterior_weights Vector containing the posterior
#'     probability of each mixture component in Ulist for the data
#'
#' @param algorithm.version Indicates whether to use R or Rcpp version
#'
#' @param A the linear transformation matrix, KxR matrix
#'
#' @return The return value is a list containing the following
#'    components:
#'
#'    \item{PosteriorMean}{J x K matrix of posterior means.}
#'
#'    \item{PosteriorSD}{J x K matrix of posterior (marginal) standard
#'    deviations.}
#'
#'    \item{NegativeProb}{J x K matrix of posterior (marginal)
#'     probability of being negative.}
#'
#'    \item{ZeroProb}{J x K matrix of posterior (marginal) probability
#'     of being zero.}
#'
#'    \item{lfsr}{J x K matrix of local false sign rates.}
#'
#' @useDynLib mashr
#'
#' @importFrom ashr compute_lfsr
#' @importFrom Rcpp evalCpp
#'
#' @export
compute_posterior_matrices <-
  function (data, Ulist, posterior_weights,
            algorithm.version = c("Rcpp","R"), A=NULL) {
  algorithm.version <- match.arg(algorithm.version)

  if(!is.null(A) && algorithm.version == 'Rcpp'){
    stop("FIXME: not implemented")
  }

  if((!is.null(data$L)) && (algorithm.version == 'Rcpp')){
    stop('FIXME: the commonbaseline method is not implemented in Rcpp')
  }

  # If A is NULL, set A be identity matrix
  R = n_conditions(data)
  if(is.null(A)){
    A = diag(R)
    row.names(A) = colnames(data$Bhat)
  }
  if(ncol(A) != R){
    stop('A is not a proper transformation')
  }

  if (algorithm.version == "R") {
    # check if covariances are same, if so, use more efficient computations
    # if alpha = 0, we check if rows of Shat are same
    # if alpha neq 0, we check if rows of Shat_alpha are same,
    # the rows of Shat_alpha are same could imply the rows of Shat are same
    common_cov_Shat = is_common_cov_Shat(data)
    if(data$alpha == 0){
      common_cov_Shat_alpha = TRUE
    } else{
      common_cov_Shat_alpha = is_common_cov_Shat_alpha(data)
    }

    if(common_cov_Shat && common_cov_Shat_alpha){ # use more efficient computations for commmon covariance case
      compute_posterior_matrices_common_cov_R(data, A, Ulist, posterior_weights)
    } else {
      compute_posterior_matrices_general_R(data, A, Ulist, posterior_weights)
    }
  } else if (algorithm.version == "Rcpp") {

    # Run the C implementation using the Rcpp interface.
    res  <- calc_post_rcpp(t(data$Bhat),t(data$Shat),data$V,
                           simplify2array(Ulist),t(posterior_weights), is_common_cov_Shat(data))
    lfsr <- compute_lfsr(res$post_neg,res$post_zero)
    posterior_matrices <- list(PosteriorMean = res$post_mean,
                              PosteriorSD   = res$post_sd,
                              lfdr          = res$post_zero,
                              NegativeProb  = res$post_neg,
                              lfsr          = lfsr)
    for (i in names(posterior_matrices)) {
      if (!is.null(colnames(data$Bhat))) colnames(posterior_matrices[[i]]) <- colnames(data$Bhat)
      if (!is.null(rownames(data$Bhat))) rownames(posterior_matrices[[i]]) <- rownames(data$Bhat)
    }
    return(posterior_matrices)
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


