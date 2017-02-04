#' Return matrix of Bhat values from mash object
#' @param m a mash object
#' @param subset.eff specifies the subset of effects (rows) to include
#' @param subset.cond specifies the subset of conditions (columns) to include
#' @export
get_Bhat = function(m, subset.eff=NULL, subset.cond=NULL){
  if(is.null(subset.eff)){subset.eff = 1:n_effects.mash(m)}
  if(is.null(subset.cond)){subset.cond = 1:n_conditions.mash(m)}
  return(m$data$Bhat[subset.eff,subset.cond])
}

#' Return matrix of Shat values from mash object
#' @param m a mash object
#' @param subset.eff specifies the subset of effects (rows) to include
#' @param subset.cond specifies the subset of conditions (columns) to include
#' @export
get_Shat = function(m, subset.eff=NULL, subset.cond=NULL){
  if(is.null(subset.eff)){subset.eff = 1:n_effects.mash(m)}
  if(is.null(subset.cond)){subset.cond = 1:n_conditions.mash(m)}
  return(m$data$Shat[subset.eff,subset.cond])
}


#' Return matrix of Z scores from mash object
#' @param m a mash object
#' @param subset.eff specifies the subset of effects (rows) to include
#' @param subset.cond specifies the subset of conditions (columns) to include
#' @export
get_Z = function(m, subset.eff=NULL, subset.cond=NULL){
  return(get_Bhat(m,subset.eff,subset.cond)/ get_Shat(m,subset.eff,subset.cond))
}



#' Return the fitted g from a mash object
#' @param m a mash object, as returned by \code{mash}
#' @export
get_fitted_g = function(m){return(list(pi=m$pi, Ulist = get_expanded_cov(m)))}

#' Return the posterior matrices from a mash object
#' @param m a mash object, as returned by \code{mash}
#' @param analysis which analysis to return results from; can be "mash" or "ash"
#' @export
get_posterior_matrices = function(m,analysis = "mash"){return(m$posterior_matrices[[analysis]])}



#' Extract grid from m
#' @param m a mash object
#' @return the grid in m
#' @export
get_grid = function(m){return(m$grid)}

#' Extract covariance matrices in m
#' @param m a mash object
#' @return a list of covariance matrices in m
#' @export
get_cov = function(m){return(m$Ulist)}

#' Get expanded list of covariance matrices in m, expanded by grid
#' @param m a mash object
#' @return a list of covariance matrices
#' This takes the covariance matrix in m and multiplies them by the grid in m
#' If a pointmass is included in m then it adds a null component
#' @export
get_expanded_cov = function(m){
  if(is.null(m$grid)){stop("need to specify grid using mash_add_grid()")}
  if(is.null(m$Ulist)){stop("need to specify some covariance matrices using add_cov()")}
  scaled_Ulist = scale_cov(m$Ulist, m$grid)
  if(m$usepointmass){scaled_Ulist = c(list(null=cov_all_zeros(m$data)),scaled_Ulist)}
  return(scaled_Ulist)
}
