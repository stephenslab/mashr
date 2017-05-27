#' Extract fitted g from mash object
#' @param m a mash object
#' @return the fitted g (a list)
#' @export
get_fitted_g = function(m){
  if(class(m)=="ash"){ashr::get_fitted_g(m)}
  return(m$fitted_g)
}

#' Extract posterior matrices from mash object
#' @param m a mash object
#' @return a list of posterior matrices computed
#' @export
get_posterior_matrices = function(m){
  return(m$result$posterior_matrices)
}

#' Extract lfdr from mash object
#' @param m a mash object
#' @return a matrix of lfsr values, with (j,r)th entry corresponding to the lfsr for effect j in condition r
#' @export
get_lfsr = function(m){
  return(ashr:::compute_lfsr(m$result$posterior_matrices$post_neg,
                             m$result$posterior_matrices$post_zero))
}

#' Find effects that have lfsr < thresh in at least one condition
#' @param m the mash result (from joint or 1by1 analysis)
#' @param thresh indicates the threshold below which to set signals
#' @return a vector containing the indices of the significant effects
#' @export
get_significant_results = function(m, thresh = 0.05){
  top_lfsr = apply(get_lfsr(m),1,min)
  which(top_lfsr< thresh)
}
