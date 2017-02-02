#' @title Perform condition-by-condition analyses
#' @param m a mash object
#' @details Performs "condition-by-condition" analysis
#' by running \code{ash} from package \code{ashr} on data from each condition, one at a time.
#' May be a useful first step to identify top hits in each condition before a mash analysis.
#' @return Adds to m the results (posterior_matrices) from the ash runs
#' @export
mash_run_1by1 = function(m, force=FALSE){
  post_mean= post_sd = lfsr = matrix(nrow = n_effects.mash(m), ncol= n_conditions.mash(m))
  for(i in 1:n_conditions.mash(m)){
    ashres = ashr::ash(get_Bhat(m,subset.cond=i),get_Shat(m,subset.cond=i),mixcompdist="normal") # get ash results for first condition
    post_mean[,i] = ashr::get_pm(ashres)
    post_sd[,i] = ashr::get_psd(ashres)
    lfsr[,i] = ashr::get_lfsr(ashres)
  }
  m$posterior_matrices[["ash"]] = list(
    post_mean = post_mean, post_sd = post_sd, lfsr = lfsr)
}
