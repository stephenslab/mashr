#' Apply mash method to data
#' @param Bhat an n by R matrix of observations (n units in R conditions)
#' @param Shat an n by R matrix of standard errors (n units in R conditions)
#' @param Ulist a list of covariance matrices to use
#' @param gridmult scalar indicating factor by which adjacent grid values should differ; close to 1 for fine grid
#' @param grid vector of grid values to use (scaling factors omega in paper)
#' @param prior indicates what penalty to use on the likelihood, if any
#' @param optmethod name of optimization method to use
#' @export
mash_new = function(Bhat,Shat,
                Ulist,
                gridmult= sqrt(2),
                grid = NULL,
                normalizeU = TRUE,
                usepointmass = TRUE,
                prior="nullbiased",
                optmethod = c("mixIP")){

  optmethod = match.arg(optmethod)
  data = set_mash_data(Bhat, Shat)
  if(missing(grid)){grid = autoselect_grid(data,gridmult)}
  if(normalizeU){Ulist = normalize_Ulist(Ulist)}
  Ulist = expand_cov(Ulist,grid,usepointmass)


  lik_matrix = calc_relative_lik_matrix(data,Ulist)$lik_matrix
  prior = set_prior(ncol(lik_matrix),prior)

  pi = optimize_pi(lik_matrix,prior=prior, optmethod=optmethod)
  posterior_weights = compute_posterior_weights(pi, lik_matrix)
  posterior_matrices = compute_posterior_matrices(data, get_expanded_cov(m), m$posterior_weights)
  return(list(posterior_matrices = posterior_matrices))

}

#' sets prior to be a vector of length K depending on character string
#' prior can be "nullbiased" or "uniform"
set_prior = function(K,prior){
  if(is.character(prior)){
    if(prior=="uniform"){
      prior=rep(1,K)
    } else if(prior=="nullbiased"){
      prior=rep(1,K); prior[1]=10
    }
  } else if(length(prior)!=K){
    stop("prior is wrong length")
  }
  return(prior)
}

#' Create expanded list of covariance matrices expanded by grid, Sigma_{lk} = omega_l U_k
#' @param Ulist a list of covarance matrices
#' @param grid a grid of scalar values by which the covariance matrices are to be sc
#' @param usepointmass if TRUE adds a point mass at 0 (null component) to the list
#' @return A list of covariance matrices
#' This takes the covariance matrices in Ulist and multiplies them by the grid values
#' If usepointmass is TRUE then it adds a null component.
#' @export
expand_cov = function(Ulist,grid,usepointmass=TRUE){
  scaled_Ulist = scale_cov(Ulist, grid)
  if(usepointmass){scaled_Ulist = c(list(null=cov_all_zeros(m$data)),scaled_Ulist)}
  return(scaled_Ulist)
}
