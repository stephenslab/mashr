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
                 cov_methods = c("identity","singletons","equal_effects","simple_het"),
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

  mash_fit_g(m,prior=prior, normalizeU)
  mash_compute_posterior(m)

  return(m)
}

#' Summarize status of a mash analysis
#' @param m a mash object
#' @details Provides a summary of what has been done so far in a mash object; may provide
#' suggestions for next step
#' @export
mash_status = function(m=NULL){
  if(class(m)!="mash"){message("m is not a mash object; SUGGEST initialize mash using mash_init"); return()}
  if(is.null(m$Ulist)){message("m has data but no covariance matrices; SUGGEST mash_add_cov"); return()}
  message("m has covariance matrices with names: ",paste(names(m$Ulist),collapse=","))
  if(is.null(m$grid)){message("m has no grid; SUGGEST add grid with mash_add_grid"); return()}
  message("m has grid: ",paste(m$grid,collapse=","))
  if(!fitted(m)){message("m has not yet had the hierarchical model fit; SUGGEST fit using mash_fit_g"); return()}

  if(is.null(m$posterior_matrices["mash"])){message("m has hierarchical model fit, but no posteriors computed;
                                                     SUGGEST use mash_compute_posterior")}
  else {message("m has had hierarchical model fit and posteriors are computed")}

  if(is.null(m$posterior_matrices[["ash"]])){message("m has no 1-by-1 (condition-by-condition) analyses run;
            SUGGEST: mash_run_1by1()"); return();}
  else {message("m has had 1-by-1 analyses run")}

  if(is.null(m$strong_signals)){message("SUGGEST: adding strong signals to allow data driven covariances")}
  if(!is.null(m$strong_signals)){message("mash has ", paste(length(m$strong_signals)), " strong signals available for estimating data-driven covariances")}


}



#' Fit the Empirical Bayes model for a mash object
#' @param m a mash object
#' @param prior The weights used in the penalized likelihood
#' @param normalize_Ulist Whether to normalized covariance matrices before fitting
#' @param optmethod Which routine to use for optimizing
#' @details Computes the likelihood matrix and fits the
#' mixture proportions for a mixture of multivariate normals.
#' The mash object contains the data, and the covariance matrices and grid to use
#' The result is stored in m$pi
#' See \code{mash_add_cov} and \code{mash_add_grid}
#' @export
mash_fit_g = function(m, prior = NULL, normalize_cov=TRUE, optmethod = c("mixIP")){
  if(normalize_cov){mash_normalize_Ulist(m)}

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
  if(!fitted(m)){stop("need to fit using mash_fit_g first")}
  m$posterior_weights = compute_posterior_weights(m$pi, m$lik_matrix)
  m$posterior_matrices[["mash"]] = compute_posterior_matrices(m$data, get_expanded_cov(m), m$posterior_weights)
}

#' Compute loglikelihood for fitted mash object
#' @param m a mash object
#' @param Bhat a matrix of Bhat values
#' @param Shat a matrix of corresponding Shat values
#' @details If no data are supplied, then computes log-likelihood on data in m (ie training data).
#' If data (Bhat and Shat) are supplied then computes log-likelihood under the model fit in m.
#' The latter might be useful for computing log-likelihood on a "test" set for example.
#' @export
mash_compute_loglik = function(m,Bhat=NULL, Shat=NULL){
  if(!fitted(m)){message("You need to fit mash using mash_fit_g first"); return();}
  if(is.null(Bhat)){return(sum(log(m$lik_matrix %*% m$pi)+m$lfactors))} else {
    if(is.null(Shat)){message("Please supply Shat"); return();}
    data = set_mash_data(Bhat,Shat)
    lm_res = calc_relative_lik_matrix(data,get_expanded_cov(m))
    return(sum(log(lm_res$lik_matrix %*% m$pi) + lm_res$lfactors))
  }
}

#' Tests whether m is fitted
fitted = function(m){return(!is.null(m$pi))}

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
  m$lik_matrix = NULL
  m$lfactors = NULL # the log factors that were removed from the lik_matrix before exponentiating

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
mash_add_grid = function(m,grid=NULL){
  if(is.null(grid)){
    m$grid = autoselect_grid(m$data,sqrt(2))
  } else {
    m$grid = grid
  }
}

#' Add a list of strong signals (to be used in data-drive covariance matrices)
#' @param m the mash object
#' @param thresh indicates the threshold below which to set signals
#' @details Adds the top univariate signals (lfsr<0.05)
#' @export
mash_add_strong_signals = function(m, thresh=0.05){
  if(is.null(m$posterior_matrices[["ash"]])){message("you need to mash_run_1by1 before adding strong signals")}
  top_lfsr = apply(m$posterior_matrices[["ash"]]$lfsr,1,min)
  m$strong_signals =which(top_lfsr< thresh)
  message(paste(length(m$strong_signals)), " signals set")
}


#' Automatically select
#' TODO: this is just a place-holder
autoselect_grid = function(data,gridmult){
  message("autoselect_grid is a place-holder\n")
  return(c(0.5,1,2))
}

#' Calculate the likelihood matrix for a mash object
#' @param m the mash object
#' @details Adds a J by P likelihood matrix to m, computed using the current data, grid and covariance matrices in m
#' @export
mash_calc_lik_matrix = function(m){
  lm_res = calc_relative_lik_matrix(m$data, get_expanded_cov(m))
  m$lik_matrix = lm_res$lik_matrix
  m$lfactors = lm_res$lfactors
}


#' Normalize the covariances in mash object
#' @param m the mash object
#' @details Normalizes the covariance matrices to have max diagonal element 1
#' @export
mash_normalize_Ulist = function(m){
  m$Ulist = normalize_Ulist(m$Ulist)
}

#' Estimate the mixture weights by maximum (penalized) likelihood
#' @param matrix_lik a matrix of likelihoods, where the (i,k)th entry is the probability of observation i given it came from component k of g
#' @param pi_init numeric vector specifying value from which to initialize optimization
#' @param prior numeric vector specifying prior to use in the penalized likelihood
#' @param optmethod a string, giving name of optimization function to use
#' @param control a list of parameters to be passed to optmethod
#' @return numeric vector specifying the optimal mixture weights
optimize_pi = function(matrix_lik, pi_init = NULL,
                       prior=NULL,
                       optmethod=c("mixEM","mixIP"),
                       control=list() ){
  optmethod = match.arg(optmethod)
  K = ncol(matrix_lik)
  if(missing(prior)){prior = rep(1,K)}
  if(missing(pi_init)){pi_init = initialize_pi(K)}
  assertthat::are_equal(length(pi_init),K)
  assertthat::are_equal(length(prior),K)

  library("ashr") # I didn't manage to get do.call to work without this
  res = do.call(optmethod, args= list(matrix_lik = matrix_lik, prior=prior, pi_init = pi_init,control=control))
  return(res$pihat)
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


