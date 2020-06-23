#' @title Return the Bayes Factor for each effect
#'
#' @param m the mash result (from joint or 1by1 analysis); must have
#' been computed using usepointmass=TRUE
#'
#' @return if m was fitted using usepointmass=TRUE then returns a
#' vector of the log10(bf) values for each effect. That is, the jth
#' element lbf[j] is log10(Pr(Bj | g=ghat-nonnull)/Pr(Bj | g = 0))
#' where ghat-nonnull is the non-null part of ghat.  Otherwise returns
#' NULL.
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data))
#' get_log10bf(m)
#'
#' @export
#'
get_log10bf = function(m) {
  if(is.null(m$null_loglik)){
    return(NULL)
  } else {
    return((m$alt_loglik - m$null_loglik)/log(10))
  }
}

#' @title Find effects that are significant in at least one condition
#'
#' @param m the mash result (from joint or 1by1 analysis)
#'
#' @param thresh indicates the threshold below which to call signals
#' significant
#'
#' @param conditions which conditions to include in check (default to all)
#'
#' @param sig_fn the significance function used to extract
#' significance from mash object; eg could be ashr::get_lfsr or
#' ashr::get_lfdr. (Small values must indicate significant.)
#'
#' @return a vector containing the indices of the significant effects,
#' by order of most significant to least
#'
#' @importFrom ashr get_lfsr
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data))
#' get_significant_results(m)
#' @export
#'
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

#' @title Count number of conditions each effect is significant in
#'
#' @param m the mash result (from joint or 1by1 analysis)
#'
#' @param thresh indicates the threshold below which to call signals
#' significant
#'
#' @param conditions which conditions to include in check (default to
#' all)
#'
#' @param sig_fn the significance function used to extract
#' significance from mash object; eg could be ashr::get_lfsr or
#' ashr::get_lfdr
#'
#' @return a vector containing the number of significant conditions
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data))
#' get_n_significant_conditions(m)
#' @export
#'
get_n_significant_conditions = function(m, thresh = 0.05, conditions = NULL,
                                        sig_fn = get_lfsr){
  if (is.null(conditions)) {
    conditions = 1:get_ncond(m)
  }
  return(apply(sig_fn(m)[,conditions,drop=FALSE]<thresh,1,sum))
}

#' @title Return the estimated mixture proportions
#'
#' @param m the mash result
#'
#' @param dimension indicates whether you want the mixture proportions for the covariances, grid, or all
#'
#' @return a named vector containing the estimated mixture proportions.
#'
#' @details If the fit was done with `usepointmass=TRUE` then the
#' first element of the returned vector will correspond to the null,
#' and the remaining elements to the non-null covariance matrices.
#' Suppose the fit was done with $K$ covariances and a grid of length
#' $L$. If `dimension=cov` then the returned vector will be of length
#' $K$ (or $K+1$ if `usepointmass=TRUE`).  If `dimension=grid` then
#' the returned vector will be of length $L$ (or $L+1$).  If
#' `dimension=all` then the returned vector will be of length $LK$ (or
#' $LK+1$). The names of the vector will be informative for which
#' combination each element corresponds to.
#'
#' @importFrom ashr get_fitted_g
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data))
#' get_estimated_pi(m)
#' @export
#'
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

#' @title Compute the proportion of (significant) signals shared by
#' magnitude in each pair of conditions, based on the poterior mean
#'
#' @param m the mash fit
#'
#' @param factor a number in [0,1] the factor within which effects are
#' considered to be shared
#'
#' @param lfsr_thresh the lfsr threshold for including an effect in
#' the assessment
#'
#' @param FUN a function to be applied to the estimated effect sizes
#' before assessing sharing. The most obvious choice beside the
#' default 'FUN=identity' would be 'FUN=abs' if you want to ignore the
#' sign of the effects when assesing sharing.
#'
#' @details For each pair of tissues, first identify the effects that
#' are significant (by lfsr<lfsr_thresh) in at least one of the two
#' tissues. Then compute what fraction of these have an estimated
#' (posterior mean) effect size within a factor `factor` of one
#' another. The results are returned as an R by R matrix.
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data))
#' get_pairwise_sharing(m) # sharing by magnitude (same sign)
#' get_pairwise_sharing(m, factor=0) # sharing by sign
#' get_pairwise_sharing(m, FUN=abs) # sharing by magnitude when sign is ignored
#' @export
get_pairwise_sharing = function(m, factor=0.5, lfsr_thresh=0.05, FUN= identity){
  R = get_ncond(m)
  lfsr = get_lfsr(m)
  S=matrix(NA,nrow = R, ncol=R)
  for(i in 1:R){
    for(j in i:R){
      sig_i=get_significant_results(m,thresh=lfsr_thresh,conditions = i)
      sig_j=get_significant_results(m,thresh=lfsr_thresh,conditions = j)
      a=union(sig_i,sig_j)
      ratio=FUN(get_pm(m)[a,i])/FUN(get_pm(m)[a,j])##divide effect sizes
      S[i,j]=mean(ratio>factor & ratio<(1/factor))
    }
  }
  S[lower.tri(S, diag = FALSE)] = t(S)[lower.tri(S, diag = FALSE)]
  colnames(S) = row.names(S) = colnames(m$result$PosteriorMean)

  return(S)
}

#' Return samples from a mash object
#'
#' @param m The mash fit.
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data), posterior_samples=5, algorithm='R')
#' get_samples(m)
#' @export
#'
get_samples = function(m){
  m$result$PosteriorSamples
}

#' @title Compute the proportion of (significant) signals shared by
#' magnitude in each pair of conditions
#'
#' @param m the mash fit with samples from posteriors
#'
#' @param factor a number in [0,1] the factor within which effects are
#' considered to be shared
#'
#' @param lfsr_thresh the lfsr threshold for including an effect in
#' the assessment
#'
#' @param FUN a function to be applied to the estimated effect sizes
#' before assessing sharing. The most obvious choice beside the
#' default 'FUN=identity' would be 'FUN=abs' if you want to ignore the
#' sign of the effects when assesing sharing.
#'
#' @details For each pair of conditions, compute the fraction of
#' effects that are within a factor `factor` of one another. The
#' results are returned as an R by R matrix.
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data), posterior_samples=5, algorithm='R')
#' get_pairwise_sharing_from_samples(m) # sharing by magnitude (same sign)
#' get_pairwise_sharing_from_samples(m, factor=0) # sharing by sign
#' get_pairwise_sharing_from_samples(m, FUN=abs) # sharing by magnitude when sign is ignored
#' @export
get_pairwise_sharing_from_samples = function(m, factor=0.5, lfsr_thresh=0.05, FUN= identity){
  samples = get_samples(m)
  if(is.null(samples)){
    stop('There is no sample from posteriors! Please use get_pairwise_sharing.')
  }
  M = dim(samples)[3]
  R = get_ncond(m)
  S = matrix(NA,nrow = R, ncol=R)
  for(i in 1:R){
    for(j in i:R){
      sig_i=get_significant_results(m,thresh=lfsr_thresh,conditions = i)
      sig_j=get_significant_results(m,thresh=lfsr_thresh,conditions = j)
      a=union(sig_i,sig_j)
      ratio = FUN(samples[a,i,])/FUN(samples[a,j,])
      S[i,j] = mean(apply(ratio, 1, function(x) sum(x > factor & x < (1/factor), na.rm = TRUE)/M))
    }
  }
  S[lower.tri(S, diag = FALSE)] = t(S)[lower.tri(S, diag = FALSE)]
  colnames(S) = row.names(S) = colnames(m$result$PosteriorMean)
  return(S)
}

get_ncond = function(m){
  return(ncol(get_pm(m)))
}