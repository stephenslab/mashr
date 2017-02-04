# notes on ED
# 1. clone the ED repo
# 2a. make rpackage (after brew install gsl ?) [not sure if necessary]
# 2b remove -fopenmp flags from r/src/Makefile
# 3. install.packages("~/Dropbox/Documents/git/extreme-deconvolution/r",repos=NULL,type="src")

#' Fit extreme deconvolution to mash data
#' @param data mash data object
#' @param Ulist_init a list of covariance matrices to initialize to
#' @return the fitted mixture: a list of mixture proportions and covariance matrices
#' @details This is a wrapper to ExtremeDeconvolution::extreme_deconvolution
#' It fixes the projection to be the identity, and the means to be 0
#' @export
mash_ed = function(data, Ulist_init){
  K = length(Ulist_init)
  R = n_conditions(data)
  pi_init = rep(1/K, K) # initial mix proportions
  ed.res = ExtremeDeconvolution:::extreme_deconvolution(data$Bhat,
                                                       data$Shat^2,
                                                       xamp = pi_init,
                                                       xmean = matrix(0,nrow=K,ncol=R),
                                                       xcovar = Ulist_init,
                                                       fixmean = TRUE)
  return(list(pi = ed.res$xamp, Ulist = ed.res$xcovar, av_loglik = ed.res$avgloglikedata))
}
