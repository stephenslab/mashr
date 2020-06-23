#' @title posterior_cov
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @details If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This
#'   function returns U1.
#'
#' @param Vinv R x R inverse covariance matrix for the likelihood
#'
#' @param U R x R prior covariance matrix
#'
#' @return R x R posterior covariance matrix
#'
#' @keywords internal
#'
posterior_cov <- function(Vinv, U){
  return(U %*% solve(Vinv %*% U + diag(nrow(U))))
}

#' @title posterior_mean
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @details If bhat is N(b,V) and b is N(0,U) then b|bhat
#'   N(mu1,U1). This function returns mu1.
#'
#' @param bhat R vector of observation
#'
#' @param Vinv R x R inverse covariance matrix for the likelihood
#'
#' @param U1 R x R posterior covariance matrix, computed using posterior_cov
#'
#' @return R vector of posterior mean
#'
#' @keywords internal
#'
posterior_mean <- function(bhat, Vinv, U1){
  return(U1 %*% (Vinv %*% bhat))
}

#' @title posterior_mean_matrix
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @details Computes posterior mean under multivariate normal model
#'   for each row of matrix Bhat. Note that if bhat is N_R(b,V) and b
#'   is N_R(0,U) then b|bhat N_R(mu1,U1). This function returns a
#'   matrix with jth row equal to mu1(bhat) for bhat= Bhat[j,].
#'
#' @param Bhat J by R matrix of observations
#'
#' @param Vinv R x R inverse covariance matrix for the likelihood
#'
#' @param U1 R x R posterior covariance matrix, computed using
#'   posterior_cov
#'
#' @return R vector of posterior mean
#'
#' @keywords internal
#'
posterior_mean_matrix <- function(Bhat, Vinv, U1){
  return(Bhat %*% (Vinv %*% U1))
}

#' @title Compute posterior matrices.
#'
#' @param data A \code{mash} data object; e.g., created by
#'     \code{\link{mash_set_data}} or \code{\link{mash_set_data_contrast}}.
#'
#' @param Ulist List containing the prior covariance matrices.
#'
#' @param posterior_weights Vector containing the posterior
#'     probability of each mixture component in Ulist for the data.
#'
#' @param algorithm.version Indicates whether to use R or Rcpp version.
#'
#' @param A the linear transformation matrix, Q x R matrix. This is
#' used to compute the posterior for Ab.
#'
#' @param mc.cores The argument supplied to
#'     \code{openmp} specifying the number of cores
#'     to use. Note that this is only has an effect for the Rcpp version.
#'
#' @param output_posterior_cov whether or not to output posterior
#' covariance matrices for all effects.
#'
#' @param posterior_samples the number of samples to be drawn from the
#' posterior distribution of each effect.
#'
#' @param seed a random number seed to use when sampling from the
#' posteriors. It is used when \code{posterior_samples > 0}.
#'
#' @return The return value is a list containing the following
#'    components:
#'
#'    \item{PosteriorMean}{J x Q matrix of posterior means.}
#'
#'    \item{PosteriorSD}{J x Q matrix of posterior (marginal) standard
#'    deviations.}
#'
#'    \item{NegativeProb}{J x Q matrix of posterior (marginal)
#'     probability of being negative.}
#'
#'    \item{ZeroProb}{J x Q matrix of posterior (marginal) probability
#'     of being zero.}
#'
#'    \item{lfsr}{J x Q matrix of local false sign rates.}
#'
#'    \item{PosteriorCov}{Q x Q x J array of posterior covariance
#'      matrices, if the \code{output_posterior_cov = TRUE}.}
#'
#'    \item{PosteriorSamples}{M x Q x J array of samples, if the
#'      \code{posterior_samples = M > 0}.}
#'
#' @useDynLib mashr
#'
#' @importFrom ashr compute_lfsr
#' @importFrom Rcpp evalCpp
#'
#' @keywords internal
#'
compute_posterior_matrices <-
  function (data, Ulist, posterior_weights,
            algorithm.version = c("Rcpp","R"), A=NULL, output_posterior_cov=FALSE,
            mc.cores = 1,
            posterior_samples = 0, seed = 123) {
  algorithm.version <- match.arg(algorithm.version)

  R = n_conditions(data)
  if (output_posterior_cov) output_type = 4
  else output_type = 3
  is_null_A = FALSE # for use with Cpp version
  # In the commonbaseline model, if the reference condition is mean, we recover the deleted column.
  if(!is.null(data$L) && attr(data$L, "reference") == 'mean'){
    temp = diag(R)
    temp = rbind(temp, -1)
    row.names(temp) = paste0(colnames(data$L), '-mean')
    if(!is.null(A)){
      if(ncol(A) != R+1){
        stop('Reference:mean. A is not a proper transformation')
      }
      A = A %*% temp
    }else{
      A = temp
    }
  }else{
    if(is.null(A)){
      A = diag(R)
      row.names(A) = colnames(data$Bhat)
      is_null_A = TRUE
    }
    if(ncol(A) != R){
      stop('A is not a proper transformation')
    }
  }
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
  is_common_cov = common_cov_Shat && common_cov_Shat_alpha
  # keep data dimension names
  effect_names = rownames(data$Bhat)
  condition_names = rownames(A)
  if (algorithm.version == "R") {
    if(data$commonV && is_common_cov){ # use more efficient computations for commmon covariance case
      posterior_matrices = compute_posterior_matrices_common_cov_R(data, A, Ulist, posterior_weights, output_posterior_cov,
                                              posterior_samples = posterior_samples, seed=seed)
    }else {
      posterior_matrices = compute_posterior_matrices_general_R(data, A, Ulist, posterior_weights, output_posterior_cov,
                                           posterior_samples = posterior_samples, seed=seed)
    }
  } else if (algorithm.version == "Rcpp") {
    if(posterior_samples > 0){
      stop('The sampling method is not implemented in C++. Please use option algorithm = "R".')
    }
    if(!data$commonV){
      stop('effect specific V has not implemented in Rcpp')
    }
    # Run the C implementation using the Rcpp interface.
    if (is_null_A) A = matrix(0,0,0)
    if (is.null(data$L))
      res <- calc_post_rcpp(t(data$Bhat), t(data$Shat), t(data$Shat_alpha), matrix(0,0,0),
                           data$V, matrix(0,0,0), A,
                           simplify2array(Ulist), t(posterior_weights),
                           is_common_cov, output_type, mc.cores)
    else
      res <- calc_post_rcpp(t(data$Bhat), t(data$Shat), t(data$Shat_alpha), t(data$Shat_orig),
                           data$V, data$L, A,
                           simplify2array(Ulist), t(posterior_weights),
                           is_common_cov, output_type, mc.cores)
    lfsr <- compute_lfsr(res$post_neg, res$post_zero)
    posterior_matrices <- list(PosteriorMean = res$post_mean,
                              PosteriorSD   = res$post_sd,
                              lfdr          = res$post_zero,
                              NegativeProb  = res$post_neg,
                              lfsr          = lfsr)

    if (output_posterior_cov) {
      posterior_matrices$PosteriorCov <- res$post_cov
    }
  } else {
    stop("Algorithm version should be either \"R\" or \"Rcpp\"")
  }
  # Set dimension names
  for (i in names(posterior_matrices)) {
    if (length(dim(posterior_matrices[[i]])) == 2) {
      colnames(posterior_matrices[[i]]) <- condition_names
      rownames(posterior_matrices[[i]]) <- effect_names
    }
  }
  if (length(dim(posterior_matrices$PosteriorCov)) == 3)
    dimnames(posterior_matrices$PosteriorCov) <- list(condition_names, condition_names, effect_names)
  return(posterior_matrices)
}

#' @title Compute posterior probabilities that each effect came from
#'   each component
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param pi a K vector of mixture proportions
#'
#' @param lik_mat a JxK matrix of likelihoods
#'
#' @return a JxK matrix of posterior probabilities, the jth row
#'   contains posteriors for jth effect
#'
#' @keywords internal
#'
compute_posterior_weights <- function(pi, lik_mat) {
  d    <- t(pi * t(lik_mat))
  norm <- rowSums(d) # normalize probabilities to sum to 1
  return(d/norm)
}
