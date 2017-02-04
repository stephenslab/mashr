#todo
#  implement possible filter of data before data-driven covs?
# implement data-driven covs (pca)

#' Apply mash method to data
#' @param Bhat an n by R matrix of observations (n units in R conditions)
#' @param Shat an n by R matrix of standard errors (n units in R conditions)
#' @param cov_methods a string indicating what covariance matrix methods to use
#' @param gridmult scalar indicating factor by which adjacent grid values should differ; close to 1 for fine grid
#' @param grid vector of grid values to use (scaling factors omega in paper)
#' @param prior indicates what penalty to use on the likelihood
#' @param optmethod name of optimization method to use
#' @export
mash = function(Bhat,Shat,
                 cov_methods = c("identity","singletons","all_ones","simple_het"),
                 gridmult= sqrt(2),
                 grid = NULL,
                 normalizeU = TRUE,
                 prior="nullbiased",
                 optmethod = c("mixIP")){

  optmethod = match.arg(optmethod)
  if(missing(grid)){grid = autoselect_grid(data,gridmult)}
  #filtered_data = filter_mash_data(data) extract top Z scores

  m=mash_init(Bhat,Shat)
  mash_add_cov(m,cov_methods) # Set up covariances

  mash_add_grid(m,grid)

  # compute likelihood matrix and optimize mixture proportions
  if(normalizeU){mash_normalize_Ulist(m)}

  mash_fit_g(m,prior=prior)
  mash_compute_posterior(m)

  return(m)
}

#' Fit the Empirical Bayes model for a mash object
#' @param m a mash object
#' @details Computes the likelihood matrix and fits the
#' mixture proportions for a mixture of multivariate normals.
#' The mash object contains the data, and the covariance matrices and grid to use
#' See \code{mash_add_cov} and \code{mash_add_grid}
#' @export
mash_fit_g = function(m, prior = NULL, optmethod = c("mixIP")){
  mash_calc_lik_matrix(m)
  K = ncol(m$lik_matrix)
  if(missing(prior)){prior="nullbiased"} #default
  if(is.character(prior)){
    if(prior=="uniform"){prior=rep(1,K)}
    else if(prior=="nullbiased"){prior=rep(1,K); prior[1]=10}
  }

  m$optdetails = list(prior = prior, optmethod = optmethod)
  m$pi=optimize_pi(m$lik_matrix, prior=prior, optmethod=optmethod)
}

#' Compute posterior quantities for mash object
#' @param m a mash object
#' @details Computes posterior weights and posterior matrices (eg posterior mean etc)
#' Uses data (Bhat and Shat) added to m and the fitted g obtained by \code{mash_fit_g}
#' @export
mash_compute_posterior = function(m){
  if(is.null(m$pi)){stop("need to fit using mash_fit_g first")}
  m$posterior_weights = compute_posterior_weights(m$pi, m$lik_matrix)
  m$posterior_matrices[["mash"]] = compute_posterior_matrices(m$data, get_expanded_cov(m), m$posterior_weights)
}

#' Initialize a mash object (actually an environment)
#' @param R the number of conditions to be used
#' @return a mash object
#' @export
mash_init = function(Bhat,Shat,usepointmass=TRUE){
  m = new.env()
  m$data = set_mash_data(Bhat, Shat)
  m$Ulist = NULL
  m$grid = NULL
  m$pi = NULL #this is currently used to check if optimized... may want to update this
  m$usepointmass = usepointmass # default is to use pointmass
  m$posterior_matrices = list()
  class(m) = "mash"
  return(m)
}

#' Add a grid of scaling factors to a mash object
#' @param m the mash object
#' @param grid a vector of numeric scaling factors
#' @details These scaling factors are used to scale all
#' covariance matrices in m when fitting the model.
#' They correspond to the omega values in Urbut et al.
#' @export
mash_add_grid = function(m,grid){
  m$grid = grid
}

#' Calculate the likelihood matrix for a mash object
#' @param m the mash object
#' @details Adds an n by R likelihood matrix to m, computed using the current data, grid and covariance matrices in m
#' @export
mash_calc_lik_matrix = function(m){
  m$lik_matrix = calc_relative_lik_matrix(m$data, get_expanded_cov(m))
}


#' Normalize the covariances in mash object
#' @param m the mash object
#' @details Normalizes the covariance matrices to have max diagonal element 1
#' @export
mash_normalize_Ulist = function(m){
  m$Ulist = normalize_Ulist(m$Ulist)
}



#' List names of covariance matrices in m
#' @param m a mash object
#' @return names of covariance matrices in m
#' @export
list_cov = function(m){names(get_cov(m))}

#' @export
n_conditions.mash = function(m){return(n_conditions(m$data))}
#' @export
n_effects.mash = function(m){return(n_effects(m$data))}


# #' Print out the components with largest weight (those exceeding thresh)
# #' @param m a mash object
# #' @param thresh the threshold on mixture weight; only components exceeding weight are output
# #' @export
# print_biggest_comp = function(m,thresh=0.01){
#   subset = which(m$pi>thresh)
#   o = order(m$pi[subset],decreasing=TRUE)
#   print(m$pi[subset][o],digits=2)
#   print(names(m$Ulist)[subset][o])
# }


