#' @title Estimate null correlations
#' @description Estimates a null correlation matrix from data
#' @param data a mash data object, eg as created by \code{mash_set_data}
#' @param Ulist a list of covariance matrices to use
#' @param init the initial value for the null correlation. If it is not given, we use result from \code{estimate_null_correlation_adhoc}
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param est_cor whether to estimate correlation matrix (TRUE) or the covariance matrix (FALSE)
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param prior indicates what penalty to use on the likelihood, if any
#' @param ... other parameters pass to \code{mash}
#' @details Returns the estimated correlation matrix (or covariance matrix) among conditions under the null.
#' The correlation (or covariance) matrix is estimated by maximum likelihood.
#' Specifically, the unknown correlation/covariance matrix V and the unknown weights are estimated iteratively.
#' The unknown correlation/covariance matrix V is estimated using M step from the EM algorithm.
#' The unknown weights pi is estimated by maximum likelihood, which is a convex problem.
#'
#' Warning: This method could take some time.
#' The \code{\link{estimate_null_correlation_simple}} gives a quick approximation for the null correlation (or covariance) matrix.
#'
#' @return a list with estimated correlation (or covariance) matrix and the fitted mash model
#' \cr
#' \item{V}{estimated correlation (or covariance) matrix}
#' \item{mash.model}{fitted mash model}
#' @importFrom stats cov2cor
#'
#' @export
#'
estimate_null_correlation_different = function(data, Ulist, init, max_iter = 30, tol=1,
                                               est_cor = TRUE, track_fit = FALSE, prior = c('nullbiased', 'uniform'),...){
  if(class(data) != 'mash'){
    stop('data is not a "mash" object')
  }
  if(!is.null(data$L)){
    stop('We cannot estimate the null correlation for the mash contrast data.')
  }

  prior <- match.arg(prior)
  tracking = list()

  if(missing(init)){
    init = tryCatch(estimate_null_correlation_simple(data, est_cor = est_cor), error = function(e) FALSE)
    if(is.logical(init)){
      warning('Use Identity matrix as the initialize null correlation.')
      init = diag(n_conditions(data))
    }
  }

  J = n_effects(data)
  m.model = fit_mash_V(data, Ulist, V = init, prior=prior,...)
  pi_s = get_estimated_pi(m.model, dimension = 'all')
  prior.v <- set_prior(length(pi_s), prior)

  # compute loglikelihood
  log_liks <- numeric(max_iter+1)
  log_liks[1] <- get_loglik(m.model) #+penalty(prior.v, pi_s)
  V = init

  result = list(V = V, mash.model = m.model)

  niter = 0
  while (niter < max_iter){
    niter = niter + 1
    if(track_fit){
      tracking[[niter]] = result
    }
    # max_V
    browser()
    V = E_V_different(data, m.model, est_cor = est_cor)

    m.model = fit_mash_V(data, Ulist, V, prior=prior, ...)
    pi_s = get_estimated_pi(m.model, dimension = 'all')

    log_liks[niter+1] <- get_loglik(m.model) # +penalty(prior.v, pi_s)

    delta.ll <- log_liks[niter+1] - log_liks[niter]

    if(delta.ll < 0){
      break
    }

    result = list(V = V, mash.model = m.model)

    if (delta.ll <= tol){
      break
    }

  }

  log_liks = log_liks[1:(niter+1)] #remove tailing NAs
  result$loglik = log_liks
  result$niter = niter + 1
  if(track_fit){
    result$trace = tracking
  }

  return(result)
}

#' @importFrom plyr aaply laply
E_V_different = function(data, m.model, est_cor){

  J = n_effects(data)
  Z = data$Bhat/data$Shat
  Shat = data$Shat * data$Shat_alpha
  post.m.shat = m.model$result$PosteriorMean / Shat

  V.array = laply( 1:J, function(j){
    temp = tcrossprod(Z[j,]) - tcrossprod(Z[j,], post.m.shat[j,]) -
      tcrossprod(post.m.shat[j,], Z[j,]) + t(m.model$result$PosteriorCov[,,j]/Shat[j,])/Shat[j,] + tcrossprod(post.m.shat[j,])
    if(est_cor){
      temp = cov2cor(temp)
    }
    return(temp)
  })

  aperm(V.array, c(2,3,1))
}

