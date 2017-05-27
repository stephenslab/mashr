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

  # calculate likelihood matrix
  lm = calc_relative_lik_matrix(data,Ulist)

  # set up prior
  prior = set_prior(ncol(lm$lik_matrix),prior)

  # main fitting procedure
  pi = optimize_pi(lm$lik_matrix,prior=prior, optmethod=optmethod)

  # compute posterior matrices
  posterior_weights = compute_posterior_weights(pi, lm$lik_matrix)
  posterior_matrices = compute_posterior_matrices(data, get_expanded_cov(m), m$posterior_weights)

  # compute log-likehood achieved
  loglik = compute_loglik_from_matrix_and_pi(pi,lm)

  return(list(posterior_matrices = posterior_matrices, loglik = loglik))

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

#' @title Perform condition-by-condition analyses
#' @param Bhat an n by R matrix of observations (n units in R conditions)
#' @param Shat an n by R matrix of standard errors (n units in R conditions)
#' @details Performs simple "condition-by-condition" analysis
#' by running \code{ash} from package \code{ashr} on data from each condition, one at a time.
#' May be a useful first step to identify top hits in each condition before a mash analysis.
#' @return posterior_matrices from the ash runs
#' @export
mash_run_1by1_new = function(Bhat,Shat){
  post_mean= post_sd = lfsr = matrix(nrow = nrow(Bhat), ncol= ncol(Bhat))
  for(i in 1:ncol(Bhat)){
    ashres = ashr::ash(Bhat[,i],Shat[,i],mixcompdist="normal") # get ash results for first condition
    post_mean[,i] = ashr::get_pm(ashres)
    post_sd[,i] = ashr::get_psd(ashres)
    lfsr[,i] = ashr::get_lfsr(ashres)
  }
  posterior_matrices = list(post_mean = post_mean, post_sd = post_sd, lfsr = lfsr)
  return(list(posterior_matrices=posterior_matrices))
}

#' Find effects that have lfsr < thresh in at least one condition
#' @param m the mash result (from joint or 1by1 analysis)
#' @param thresh indicates the threshold below which to set signals
#' @return a vector containing the indices of the significant effects
#' @export
get_significant_results = function(m, thresh = 0.05){
  top_lfsr = apply(m$posterior_matrices$lfsr,1,min)
  which(top_lfsr< thresh)
}

#' Compute loglikelihood from a matrix of log-likelihoods and fitted pi
#' @param pi the vector of mixture proportions
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' @export
compute_loglik_from_matrix_and_pi = function(pi,lm){
  return(sum(log(lm$lik_matrix %*% pi)+lm$lfactors))
}
