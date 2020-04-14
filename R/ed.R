# notes on ED
# 1. clone the ED repo
# 2 remove -fopenmp flags from r/src/Makefile
# 3. make rpackage (after brew install gsl ?) [not sure if necessary]
# 3alt. or possibly install.packages("~/Dropbox/Documents/git/extreme-deconvolution/r",repos=NULL,type="src")

#' Fit extreme deconvolution to mash data using Bovy et al 2011
#' @param data mash data object
#' @param Ulist_init a list of covariance matrices to initialize to
#' @param subset the indices of the observations to be used (defaults to all of them)
#' @param ... arguments to be passed to \code{extreme_deconvolution} function, such as \code{tol}, \code{maxiter}.
#' @return the fitted mixture: a list of mixture proportions and covariance matrices
#' @details This is a wrapper to ExtremeDeconvolution::extreme_deconvolution
#' It fixes the projection to be the identity, and the means to be 0
#'
#' @export
bovy_wrapper = function(data, Ulist_init, subset=NULL, ...){
  if(is.null(subset)){subset = 1:n_effects(data)}
  K = length(Ulist_init)
  R = n_conditions(data)
  pi_init = rep(1/K, K) # initial mix proportions
  D = ncol(data$V)
  if(all(data$V==diag(D))){
    ed.res = extreme_deconvolution(data$Bhat[subset,],
                                   data$Shat[subset,]^2,
                                   xamp = pi_init,
                                   xmean = matrix(0,nrow=K,ncol=R),
                                   xcovar = Ulist_init,
                                   fixmean = TRUE,
                                   ...)
  }else{
    ycovar = lapply(subset, function(i) data$Shat[i,] * t(data$V * data$Shat[i,]) )
    ed.res = extreme_deconvolution(data$Bhat[subset,],
                                   ycovar,
                                   xamp = pi_init,
                                   xmean = matrix(0,nrow=K,ncol=R),
                                   xcovar = Ulist_init,
                                   fixmean = TRUE,
                                   ...)
  }
  return(list(pi = ed.res$xamp, Ulist = ed.res$xcovar, av_loglik = ed.res$avgloglikedata))
}

#' Fit extreme deconvolution to mash data using TEEM method developed by Y. Yang and M Stephens
#' @param data mash data object
#' @param Ulist_init a list of covariance matrices to initialize to
#' @param subset the indices of the observations to be used (defaults to all of them)
#' @return the fitted mixture: a list of mixture proportions and covariance matrices
#'
#' @export
teem_wrapper = function(data, Ulist_init, subset=NULL, w_init=NULL, maxiter=5000, tol=1e-7, verbose=FALSE) {
  if(is.null(subset)){subset = 1:n_effects(data)}
  zscore = data$Bhat[subset,]/data$Shat[subset,]
  if(is.null(w_init)) w_init = rep(1/length(Ulist_init), length(Ulist_init))
  for (i in 1:length(Ulist_init)) Ulist_init[[i]] = Ulist_init[[i]] + diag(nrow(Ulist_init[[i]]))
  res = fit_teem_rcpp(zscore, w_init, simplify2array(Ulist_init), maxiter, tol, verbose)
  # format result to list
  res$U = lapply(seq(dim(res$U)[3]), function(x) res$U[ , , x])
  return(res)
}