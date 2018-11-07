#' @title Estimate null correlations (ad hoc)
#' @description Estimates a null correlation matrix from data using simple z score threshold
#' @param data a mash data object, eg as created by \code{mash_set_data}
#' @param z_thresh the z score threshold below which to call an effect null
#' @details Returns the empirical correlation matrix of the effects that are "null" based on simple z score threshold
#' @importFrom stats cor
#' @export
estimate_null_correlation_adhoc = function(data, z_thresh=2){
  z = data$Bhat/data$Shat
  max_absz = apply(abs(z),1, max)
  nullish = which(max_absz < z_thresh)
  if(length(nullish)<n_conditions(data)){
    stop("not enough null data to estimate null correlation")
  }
  nullish_z = z[nullish,]
  Vhat = cor(nullish_z)
  return(Vhat)
}

# #' @title compute_null_correlation_lower_bound
# #' @description Compute a lower bound on the null correlation matrix from data
# #' @details Assume $$z^j_r = mu^j_r + e_r$$ where $e_1,e_2$ are joint normal with
# #' variance 1 and some covariance (same as correlation since they are variance 1).
# #' Then $E((z^j_1 - z^j_2)^2) = E(mu^j_1-mu^j_2)^2 + 2(1-cov(e_1,e_2)) > 2(1-cov(e_1,e_2))$.
# #' Thus $$cov(e_1,e_2) > 1- 0.5E((z^j_1 - z^j_2)^2)$$ gives a lower bound on the covariance.
# #'
# #' @param data Description of this argument goes here.
# compute_null_correlation_lower_bound = function(data){
#   R = n_conditions(data)
#   z = data$Bhat/data$Shat
#   lb = matrix(0,ncol=R,nrow=R) # the lower bound
#   for(i in 1:(R-1)){
#     for(j in (i+1):R){
#       lb[i,j] = 1-0.5*mean((z[,i]-z[,j])^2)
#       lb[j,i] = lb[i,j] #lb is symmetric
#     }
#   }
#   return(lb)
# }

#' @title Estimate null correlations
#' @description Estimates a null correlation matrix from data
#' @param data a mash data object, eg as created by \code{mash_set_data}
#' @param Ulist a list of covariance matrices to use
#' @param init the initial value for the null correlation. If it is not given, we use result from \code{estimate_null_correlation_adhoc}
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param cor whether to estimate correlation matrix. If it is False, we estimate the covairance matrix
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param prior indicates what penalty to use on the likelihood, if any
#' @param ... other parameters pass to \code{mash}
#' @details Returns the estimated correlation/covariance matrix of the effects
#'
#' @importFrom stats cov2cor
#' 
#' @export
#'
estimate_null_correlation = function(data, Ulist, init, max_iter = 50, tol=1e-3,
                                     cor = TRUE, track_fit = FALSE, prior = c('nullbiased', 'uniform'),...){
  if(class(data) != 'mash'){
    stop('data is not a "mash" object')
  }
  if(!is.null(data$L)){
    stop('We cannot estimate the null correlation for the mash contrast data.')
  }

  prior <- match.arg(prior)
  tracking = list()

  if(missing(init)){
    init = tryCatch(estimate_null_correlation_adhoc(data), error = function(e) FALSE)
    if(is.logical(init)){
      warning('Use Identity matrix as the initialize null correlation.')
      init = diag(n_conditions(data))
    }
  }

  m.model = fit_mash_V(data, Ulist, V = init, prior=prior)
  pi_s = get_estimated_pi(m.model, dimension = 'all')
  prior.v <- set_prior(length(pi_s), prior)

  # compute loglikelihood
  log_liks <- numeric(max_iter+1)
  log_liks[1] <- get_loglik(m.model)+penalty(prior.v, pi_s)
  V = init

  result = list(V = V, mash.model = m.model)

  niter = 0
  while (niter < max_iter){
    niter = niter + 1
    if(track_fit){
      tracking[[niter]] = result
    }
    # max_V
    V = E_V(data, m.model)
    if(cor){
      V = cov2cor(V)
    }
    m.model = fit_mash_V(data, Ulist, V, prior=prior)
    pi_s = get_estimated_pi(m.model, dimension = 'all')

    log_liks[niter+1] <- get_loglik(m.model)+penalty(prior.v, pi_s)

    result = list(V = V, mash.model = m.model)

    # Update delta
    delta.ll <- log_liks[niter+1] - log_liks[niter]
    if(delta.ll<=tol) break;
  }

  if(track_fit){
    result$trace = tracking
  }

  return(result)
}

penalty <- function(prior, pi_s){
  subset <- (prior != 1.0)
  sum((prior-1)[subset]*log(pi_s[subset]))
}

#' @importFrom plyr aaply laply
E_V = function(data, m.model){
  n = n_effects(data)
  Z = data$Bhat/data$Shat
  post.m.shat = m.model$result$PosteriorMean / data$Shat
  post.sec.shat = laply(1:n, function(i) (t(m.model$result$PosteriorCov[,,i]/data$Shat[i,])/data$Shat[i,]) + tcrossprod(post.m.shat[i,])) # nx2x2 array
  temp1 = crossprod(Z)
  temp2 = crossprod(post.m.shat, Z) + crossprod(Z, post.m.shat)
  temp3 = unname(aaply(post.sec.shat, c(2,3), sum))

  (temp1 - temp2 + temp3)/n
}

fit_mash_V <- function(data, Ulist, V, prior=c('nullbiased', 'uniform')){
  data$V = V
  m.model = mash(data, Ulist, prior=prior, verbose = FALSE, outputlevel = 3)
  return(m.model)
}
