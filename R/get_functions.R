#note that many of the get_ functions used in mashr (get_lfsr, get_pm etc) are defined in the ashr package

#' Find effects that have lfsr < thresh in at least one condition
#' @param m the mash result (from joint or 1by1 analysis)
#' @param thresh indicates the threshold below which to set signals
#' @param sig_fn the significance function used to extract significance from mash object; eg could be ashr::get_lfsr or ashr::get_lfdr
#' @return a vector containing the indices of the significant effects, by order of most significant to least
#' @export
get_significant_results = function(m, thresh = 0.05, sig_fn=ashr::get_lfsr){
  top = apply(sig_fn(m),1,min) #find top effect in each condition
  sig = which(top < thresh)
  ord = order(top[sig],decreasing=FALSE)
  sig[ord]
}


#' Return the estimated mixture proportions
#' @param m the mash result
#' @param dimension indicates whether you want the mixture proportions for the covariances, grid, or all
#' @return a named vector containing the estimated mixture proportions.
#' @details If the fit was done with `usepointmass=TRUE` then the first element of the returned vector will correspond to the null, and the remaining elements
#' to the non-null covariance matrices.
#' Suppose the fit was done with $K$ covariances and a grid of length $L$.
#' If `dimension=cov` then the returned vector will be of length $K$ (or $K+1$ if
#' `usepointmass=TRUE`).  If `dimension=grid` then the returned vector will be of length $L$ (or $L+1$).
#' If `dimension=all` then the returned vector will be of length $LK$ (or $LK+1$).
#' The names of the vector will be informative for which combination each element corresponds to.
#' @export
get_estimated_pi = function(m, dimension = c("cov","grid","all")){
  dimension = match.arg(dimension)
  if(dimension=="all"){
    get_estimated_pi_no_collapse(m)
  } else {
    g = get_fitted_g(m)
    pihat = g$pi
    pihat_names = NULL
    pi_null = NULL

    if(g$usepointmass){
      pihat_names=c("null",pihat_names)
      pi_null = pihat[1]
      pihat = pihat[-1]
    }

    pihat = matrix(pihat,nrow=length(g$Ulist))
    if(dimension=="cov"){
      pihat = rowSums(pihat)
      pihat_names = c(pihat_names,names(g$Ulist))
    } else if(dimension=="grid"){
      pihat = colSums(pihat)
      pihat_names = c(pihat_names,1:length(g$grid))
    }

    pihat = c(pi_null,pihat)
    names(pihat) = pihat_names
    return(pihat)
  }
}


get_estimated_pi_no_collapse = function(m){
  g = get_fitted_g(m)
  pihat = g$pi
  names(pihat) = names(expand_cov(g$Ulist, g$grid, g$usepointmass))
  pihat
}
