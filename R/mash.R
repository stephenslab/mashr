# a draft of what the interface might look like
#todo
# test posterior calculations; compare posterior means with ash estimates, especially in independent case
#  implement possible filter of data before data-driven covs?
# implement data-driven covs (pca)

#' Apply mash method to data
#' @param Bhat an n by R matrix of observations (n units in R conditions)
#' @param Shat an n by R matrix of standard errors (n units in R conditions)
#' @param cov_methods a string indicating what covariance matrix methods to use
#' @param gridmult scalar indicating factor by which adjacent grid values should differ; close to 1 for fine grid
#' @param grid vector of grid values to use (scaling factors omega in paper)
#' @param prior indicates what penalty to use on the likelihood
#' @export
mash = function(Bhat,Shat, cov_methods = c("identity","singletons"), gridmult= sqrt(2), grid = NULL, prior="nullbiased"){
  data = set_mash_data(Bhat,Shat)
  if(missing(grid)){grid = autoselect_grid(data,gridmult)}
  #filtered_data = filter_mash_data(data) extract top Z scores

  g = add_to_g(data,
               cov_methods,
               grid)
  g = add_to_g(data, "null", 1, g) # add null to g

  lik_matrix = calc_relative_lik_matrix(data, g$Ulist)
  g_opt=optimize_g(g, lik_matrix, prior="nullbiased", optmethod="mixIP")

  posterior_weights = compute_posterior_weights(get_mixprob(g_opt), lik_matrix)
  posterior_matrices = compute_posterior_matrices(data, g_opt, posterior_weights)

  return(list(data=data, fitted_g= g_opt, result=posterior_matrices))
}

#' Return the fitted g from a mash object
#' @param m a mash object, as returned by \code{mash}
#' @export
get_fitted_g = function(m){return(m$fitted_g)}

#' Return the posterior matrices from a mash object
#' @param m a mash object, as returned by \code{mash}
#' @export
get_posterior_matrices = function(m){return(m$result)}
