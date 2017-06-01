#note that many of the get_ functions used in mashr (get_lfsr, get_pm etc) are defined in the ashr package

#' Find effects that have lfsr < thresh in at least one condition
#' @param m the mash result (from joint or 1by1 analysis)
#' @param thresh indicates the threshold below which to set signals
#' @param conditions which conditions to include in check (default to all)
#' @param sig_fn the significance function used to extract significance from mash object; eg could be ashr::get_lfsr or ashr::get_lfdr
#' @return a vector containing the indices of the significant effects, by order of most significant to least
#' @importFrom ashr get_lfsr
#' @export
get_significant_results = function(m, thresh = 0.05, conditions = NULL,
    sig_fn = get_lfsr) {
  if (is.null(conditions)) {
    conditions = 1:get_ncond(m)
  }
  top = apply(sig_fn(m)[,conditions,drop=FALSE],1,min) # find top effect in each condition
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
#' @importFrom ashr get_fitted_g
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

#' Compute the proportion of (significant) signals shared by magnitude in each pair of conditions
#' @param m the mash fit
#' @param factor a number in [0,1] the factor within which effects are considered to be shared
#' @param lfsr_thresh the lfsr threshold for including an effect in the assessment
#' @param FUN a function to be applied to the estimated effect sizes before assessing sharing. The most obvious choice beside the default
#' 'FUN=identity' would be 'FUN=abs' if you want to ignore the sign of the effects when assesing sharing.
#' @details For each pair of tissues, first identify the effects that are significant (by lfsr<lfsr_thresh)
#'  in at least one of the two tissues. Then compute what fraction of these have an estimated (posterior mean) effect size within
#'  a factor `factor` of one another. The results are returned as an R by R matrix.
#' @examples
#' get_pairwise_sharing(m) # sharing by magnitude (same sign)
#' get_pairwise_sharing(m, factor=0) # sharing by sign
#' get_pairwise_sharing(m, FUN=abs) # sharing by magnitude when sign is ignored
#' @export
get_pairwise_sharing = function(m, factor=0.5, lfsr_thresh=0.05, FUN= identity){
  R = get_ncond(m)
  lfsr = get_lfsr(m)
  S=matrix(NA,nrow = R, ncol=R)
  #colnames(S)=rownames(S)=colnames(maxz)
  for(i in 1:R){
    for(j in 1:R){
      sig_i=get_significant_results(m,thresh=lfsr_thresh,conditions = i)
      sig_j=get_significant_results(m,thresh=lfsr_thresh,conditions = j)
      a=union(sig_i,sig_j)
      ratio=FUN(get_pm(m)[a,i])/FUN(get_pm(m)[a,j])##divide effect sizes
      S[i,j]=mean(ratio>factor & ratio<(1/factor))
    }
  }
  return(S)
}

get_ncond = function(m){
  return(ncol(get_pm(m)))
}
