#' Initialize a mixture of multivariate normals
#' @param data a mash data object, e.g. as created by \code{set_mash_data}
#' @param cov_methods methods to use to create covariance matrices; passed to \code{cov_methods}
#' @param grid a vector of scalars that are to be used to scale the covariance matrices
#' @param g optionally a previously-inititalized g; in this case the new g adds new covariances to this g
#' @return a list with two elements, the covariance matrices and their mixture proportions
#' The covariance matrices are

#' @export
initialize_g=function(data, cov_methods, grid, g=NULL){
  Ulist = compute_cov(data,cov_methods,g$Ulist)
  Ulist = normalize_Ulist(Ulist)
  Ulist = scale_cov(Ulist, grid)

  pi = initialize_pi(length(Ulist))
  list(Ulist=Ulist, pi=pi)
}

#' Initialize mixture proportions - currently by making them all equal
#' @param K the number of components
#' @return a vector of length K whose elements are positive and sums to 1
initialize_pi = function(K){
  return(rep(1/K,K))
}

#' Automatically select
#' TODO: this is just a place-holder
autoselect_grid = function(data){
  message("autoselect_grid is a place-holder\n")
  return(c(0.5,1,2))
}

#' Return the covariances in g
#' @param g a mixture of multivariate normals, as created for example by \code{initialize_g}
#' @export
get_cov = function(g){return(g$Ulist)}

#' Return the mixture proportions in g
#' @param g a mixture of multivariate normals, as created for example by \code{initialize_g}
#' @export
get_mixprob = function(g){return(g$pi)}

n_comp = function(g){return(length(g$pi))}

null_comp = function(g){
  which(grepl("all_zeros",names(g$Ulist)))
}

#' Print out the components of g with largest weight (those exceeding thresh)
#' @param g a mixture of multivariate normals, e.g. as created by \code{initialize_g}
#' @param thresh the threshold on mixture weight; only components exceeding weight are output
#' @export
print_biggest_comp = function(g,thresh=0.01){
  subset = which(g$pi>thresh)
  o = order(g$pi[subset],decreasing=TRUE)
  print(g$pi[subset][o],digits=2)
  print(names(g$Ulist)[subset][o])
}

#' Fit the hierarchical model by estimating the mixture weights
#' @param data a mash data object, e.g. as created by \code{set_mash_data}
#' @param g_init a mixture of multivariate normals, e.g. as created by \code{initialize_g}
#' @param prior a string saying what kind of prior to use in the penalized likelihood
#' @export
optimize_g = function(data,
                      g_init,
                      prior=c("nullbiased"),
                      optmethod=c("mixEM","mixIP"),
                      control=list() ){
  optmethod = match.arg(optmethod)
  library("ashr") # I didn't manage to get do.call to work without this
  prior= ashr:::setprior(prior, n_comp(g_init), 10, null_comp(g_init))
  matrix_lik = calc_relative_lik_matrix(data, get_cov(g_init))
  pi_init = get_mixprob(g_init)
  res = do.call(optmethod, args= list(matrix_lik = matrix_lik, prior=prior, pi_init = pi_init,control=control))
  g_init$pi = res$pihat
  return(g_init)
}
