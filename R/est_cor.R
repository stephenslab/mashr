#' @title Estimate null correlations (simple)
#'
#' @description Estimates a null correlation matrix from data using
#' simple z score threshold
#'
#' @param data a mash data object, eg as created by \code{mash_set_data}
#'
#' @param z_thresh the z score threshold below which to call an effect null
#'
#' @param est_cor whether to estimate correlation matrix (TRUE) or the
#' covariance matrix (FALSE).
#'
#' @details Returns a simple estimate of the correlation matrix (or
#' covariance matrix) among conditions under the null.  Specifically,
#' the simple estimate is the empirical correlation (or covariance)
#' matrix of the z scores for those effects that have (absolute) z
#' score < z_thresh in all conditions.
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' estimate_null_correlation_simple(data)
#'
#' @importFrom stats cov2cor
#' @importFrom stats cov
#'
#' @export
#'
estimate_null_correlation_simple = function(data, z_thresh=2, est_cor = TRUE){
  z = data$Bhat/data$Shat
  max_absz = apply(abs(z),1, max)
  nullish = which(max_absz < z_thresh)
  if(length(nullish)<n_conditions(data)){
    stop("not enough null data to estimate null correlation")
  }
  nullish_z = z[nullish,]
  Vhat = cov(nullish_z)
  if(est_cor){
    Vhat = cov2cor(Vhat)
  }
  return(Vhat)
}

#' @title Estimate null correlations
#'
#' @description Estimates a null correlation matrix from data
#'
#' @param data a mash data object, eg as created by \code{mash_set_data}
#'
#' @param Ulist a list of covariance matrices to use
#'
#' @param init the initial value for the null correlation. If it is
#' not given, we use result from
#' \code{estimate_null_correlation_adhoc}
#'
#' @param max_iter maximum number of iterations to perform
#'
#' @param tol convergence tolerance
#'
#' @param est_cor whether to estimate correlation matrix (TRUE) or the
#' covariance matrix (FALSE)
#'
#' @param track_fit add an attribute \code{trace} to output that saves
#' current values of all iterations
#'
#' @param prior indicates what penalty to use on the likelihood, if any
#'
#' @param details whether to return details of the model, if it is
#' TRUE, the number of iterations and the value of objective functions
#' will be returned
#'
#' @param ... other parameters pass to \code{mash}
#'
#' @details Returns the estimated correlation matrix (or covariance
#' matrix) among conditions under the null.  The correlation (or
#' covariance) matrix is estimated by maximum likelihood.
#' Specifically, the unknown correlation/covariance matrix V and the
#' unknown weights are estimated iteratively.  The unknown
#' correlation/covariance matrix V is estimated using the posterior
#' second moment of the noise.  The unknown weights pi is estimated by
#' maximum likelihood, which is a convex problem.
#'
#' Warning: This method could take some time.  The
#' \code{\link{estimate_null_correlation_simple}} gives a quick
#' approximation for the null correlation (or covariance) matrix.
#'
#' @return the estimated correlation (or covariance) matrix and the
#' fitted mash model \cr
#'
#' \item{V}{estimated correlation (or covariance) matrix}
#'
#' \item{mash.model}{fitted mash model}
#'
#' @examples
#' simdata = simple_sims(100,5,1)
#' m.1by1 = mash_1by1(mash_set_data(simdata$Bhat,simdata$Shat))
#' strong.subset = get_significant_results(m.1by1,0.05)
#' random.subset = sample(1:nrow(simdata$Bhat),20)
#' data.strong = mash_set_data(simdata$Bhat[strong.subset,], simdata$Shat[strong.subset,])
#' data.tmp = mash_set_data(simdata$Bhat[random.subset,], simdata$Shat[random.subset,])
#' U_pca = cov_pca(data.strong, 3)
#' U_ed = cov_ed(data.strong, U_pca)
#' Vhat = estimate_null_correlation(data.tmp, U_ed)
#' @importFrom stats cov2cor
#'
#' @export
#'
estimate_null_correlation = function(data, Ulist, init, max_iter = 30, tol=1,
                                     est_cor = TRUE, track_fit = FALSE, prior = c('nullbiased', 'uniform'), details = FALSE, ...){
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
    V = E_V(data, m.model)/J
    if(est_cor){
      V = cov2cor(V)
    }
    m.model = fit_mash_V(data, Ulist, V, prior=prior, ...)
    pi_s = get_estimated_pi(m.model, dimension = 'all')

    log_liks[niter+1] <- get_loglik(m.model)  #+penalty(prior.v, pi_s)
    delta.ll <- log_liks[niter+1] - log_liks[niter]
    if(delta.ll < 0){
      break
    }

    result = list(V = V, mash.model = m.model)
    if (delta.ll <= tol){
      niter = niter + 1
      break
    }
  }

  log_liks = log_liks[1:niter] #remove tailing NAs
  result$loglik = log_liks
  result$niter = niter
  if(track_fit){
    result$trace = tracking
  }

  if(details){
    return(result)
  }else{
    return(result$V)
  }
}

#' @importFrom plyr aaply laply
E_V = function(data, m.model){
  J = n_effects(data)
  Z = data$Bhat/data$Shat
  Shat = data$Shat * data$Shat_alpha
  post.m.shat = m.model$result$PosteriorMean / Shat
  post.sec.shat = laply(1:J, function(i) (t(m.model$result$PosteriorCov[,,i]/Shat[i,])/Shat[i,]) +
                          tcrossprod(post.m.shat[i,])) # JxRxR array
  temp1 = crossprod(Z)
  temp2 = crossprod(post.m.shat, Z) + crossprod(Z, post.m.shat)
  temp3 = unname(aaply(post.sec.shat, c(2,3), sum))

  V = (temp1 - temp2 + temp3)
  # avoid numerical unsymmetry
  V = (V+t(V))/2
}

fit_mash_V <- function(data, Ulist, V, prior=c('nullbiased', 'uniform'), ...){
  data.V = mash_update_data(data, V=V)
  m.model = mash(data.V, Ulist, prior=prior, verbose = FALSE, outputlevel = 3, ...)
  return(m.model)
}
