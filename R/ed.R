#' @title Fit extreme deconvolution to mash data using Bovy et al 2011
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param data mash data object
#'
#' @param Ulist_init a list of covariance matrices to initialize to
#'
#' @param subset the indices of the observations to be used (defaults
#' to all of them)
#'
#' @param ... arguments to be passed to \code{extreme_deconvolution}
#' function, such as \code{tol}, \code{maxiter}.
#'
#' @return the fitted mixture: a list of mixture proportions and
#' covariance matrices
#'
#' @details This is a wrapper to ExtremeDeconvolution::extreme_deconvolution
#' It fixes the projection to be the identity, and the means to be 0
#'
#' @keywords internal
#'
bovy_wrapper = function(data, Ulist_init, subset=NULL, ...){
  if(is.null(subset)){subset = 1:n_effects(data)}
  K = length(Ulist_init)
  R = n_conditions(data)
  pi_init = rep(1/K, K) # initial mix proportions
  D = ncol(data$V)
  if(!is.null(data$L)){
    ycovar = lapply(subset, function(i) data$L %*% (data$Shat_orig[i,] * t(data$V * data$Shat_orig[i,])) %*% t(data$L) )
  }else if(!all(data$V==diag(D))){
    ycovar = lapply(subset, function(i) data$Shat[i,] * t(data$V * data$Shat[i,]) )
  }else{
    ycovar = data$Shat[subset,]^2
  }
  ed.res = extreme_deconvolution(data$Bhat[subset,],
                                 ycovar,
                                 xamp = pi_init,
                                 xmean = matrix(0,nrow=K,ncol=R),
                                 xcovar = Ulist_init,
                                 fixmean = TRUE,
                                 ...)
  return(list(pi = ed.res$xamp, Ulist = ed.res$xcovar, av_loglik = ed.res$avgloglikedata))
}

#' @title Fit extreme deconvolution to mash data using TEEM method
#' developed by Y. Yang and M Stephens
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param data mash data object
#'
#' @param Ulist_init a list of covariance matrices to initialize to
#'
#' @param subset the indices of the observations to be used (defaults
#' to all of them)
#'
#' @return the fitted mixture: a list of mixture proportions and
#' covariance matrices
#'
#' @keywords internal
#'
teem_wrapper = function(data, Ulist_init, subset=NULL, w_init=NULL, maxiter=5000, converge_tol=1e-7, eigen_tol = 1e-7, verbose=FALSE) {
  if(is.null(subset)){subset = 1:n_effects(data)}
  zscore = data$Bhat[subset,]/data$Shat[subset,]
  if(is.null(w_init)) w_init = rep(1/length(Ulist_init), length(Ulist_init))
  res = fit_teem_rcpp(zscore, w_init, simplify2array(Ulist_init), maxiter, converge_tol, eigen_tol, verbose)
  # format result to list with names
  names(res$U) = names(Ulist_init)
  res$w = as.vector(res$w)
  res$objective = as.vector(res$objective)
  res$maxd = as.vector(res$maxd)
  return(res)
}
