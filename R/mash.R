#' Apply mash method to data
#' @param data a mash data object containing the Bhat matrix and standard errors
#' @param Ulist a list of covariance matrices to use
#' @param gridmult scalar indicating factor by which adjacent grid values should differ; close to 1 for fine grid
#' @param grid vector of grid values to use (scaling factors omega in paper)
#' @param normalizeU whether or not to normalize the U covariances to have maximum of 1 on diagonal
#' @param usepointmass whether to include a point mass at 0, corresponding to null in every condition
#' @param g the value of g obtained from a previous mash fit - an alternative to supplying Ulist, grid and usepointmass
#' @param fixg if g is supplied, allows the mixture proportions to be fixed rather than estimated - e.g. useful for fitting mash to test data after fitting it to training data
#' @param prior indicates what penalty to use on the likelihood, if any
#' @param optmethod name of optimization method to use
#' @return a list with elements result, loglik and fitted_g
#' @export
mash = function(data,
                Ulist = NULL,
                gridmult= sqrt(2),
                grid = NULL,
                normalizeU = TRUE,
                usepointmass = TRUE,
                g = NULL,
                fixg = FALSE,
                prior=c("nullbiased","uniform"),
                optmethod = c("mixIP","mixEM","cxxMixSquarem")){



  if(!missing(g)){ # g is supplied
    if(!missing(Ulist)){stop("cannot supply both g and Ulist")}
    if(!missing(grid)){stop("cannot supply both g and grid")}
    if(!missing(normalizeU)){stop("cannot supply both g and normalizeU")}
    if(!missing(usepointmass)){stop("cannot supply both g and usepointmass")}
    Ulist = g$Ulist
    grid = g$grid
  } else { #g not supplied
    if(missing(Ulist)){stop("must supply Ulist (or g from previous mash fit)")}
    if(missing(grid)){grid = autoselect_grid(data,gridmult)}
    if(normalizeU){Ulist = normalize_Ulist(Ulist)}
  }

  if(fixg){
    if(missing(g)){stop("cannot fix g if g not supplied!")}
    if(!missing(prior)){stop("cannot supply prior if fixg is TRUE")}
    if(!missing(optmethod)){stop("cannot supply optmethod if fixg is TRUE")}
  } else {
    optmethod = match.arg(optmethod)
    prior = match.arg(prior)
  }


  xUlist = expand_cov(Ulist,grid,usepointmass)

  # calculate likelihood matrix
  lm = calc_relative_lik_matrix(data,xUlist)

  # main fitting procedure
  if(!fixg){
    prior = set_prior(ncol(lm$lik_matrix),prior)
    pi = optimize_pi(lm$lik_matrix,prior=prior, optmethod=optmethod)
  }
  else{ #if fixg, just use g$pi for pi
    pi = g$pi
  }

  # compute posterior matrices
  posterior_weights = compute_posterior_weights(pi, lm$lik_matrix)
  posterior_matrices = compute_posterior_matrices(data, xUlist, posterior_weights)

  # compute log-likehood achieved
  loglik = compute_loglik_from_matrix_and_pi(pi,lm)
  fitted_g = list(pi = pi, Ulist=Ulist, grid=grid, usepointmass=usepointmass)

  m=list(result=posterior_matrices, loglik = loglik, fitted_g = fitted_g)
  class(m) = "mash"
  return(m)
}

#' Compute loglikelihood for fitted mash object on new data
#' @param g a mash object or the fitted_g from a mash object
#' @param data a set of data on which to compute the loglikelihood
#' @return the log-likelihood for data computed using g
#' @export
mash_compute_loglik = function(g,data){
  if(class(g)=="mash"){g = g$fitted_g}
  xUlist = expand_cov(g$Ulist,g$grid,g$usepointmass)
  lm_res = calc_relative_lik_matrix(data,xUlist)
  return(sum(log(lm_res$lik_matrix %*% g$pi) + lm_res$lfactors))
}

#' Compute posterior matrices for fitted mash object on new data
#' @param g a mash object or the fitted_g from a mash object
#' @param data a set of data on which to compute the posterior matrices
#' @return A list of posterior matrices
#' @export
mash_compute_posterior_matrices = function(g,data){
  if(class(g)=="mash"){g = g$fitted_g}

  xUlist = expand_cov(g$Ulist,g$grid,g$usepointmass)
  lm_res = calc_relative_lik_matrix(data,xUlist)

  posterior_weights = compute_posterior_weights(g$pi, lm_res$lik_matrix)
  posterior_matrices = compute_posterior_matrices(data, xUlist, posterior_weights)

  return(posterior_matrices)
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
  R = nrow(Ulist[[1]])
  if(usepointmass){
    scaled_Ulist = c(list(null=matrix(0,nrow=R,ncol=R)),scaled_Ulist)
  }
  return(scaled_Ulist)
}

#' @title Perform condition-by-condition analyses
#' @param Bhat an n by R matrix of observations (n units in R conditions)
#' @param Shat an n by R matrix of standard errors (n units in R conditions)
#' @description Performs simple "condition-by-condition" analysis
#' by running \code{ash} from package \code{ashr} on data from each condition, one at a time.
#' May be a useful first step to identify top hits in each condition before a mash analysis.
#' @return a list similar to the output of mash, particularly including posterior matrices
#' @importFrom ashr ash get_pm get_psd get_lfsr get_loglik
#' @export
mash_1by1 = function(data){
  Bhat = data$Bhat
  Shat = data$Shat
  post_mean= post_sd = lfsr = matrix(nrow = nrow(Bhat), ncol= ncol(Bhat))
  loglik = 0
  for(i in 1:ncol(Bhat)){
    ashres = ash(Bhat[,i],Shat[,i],mixcompdist="normal") # get ash results for first condition
    post_mean[,i] = get_pm(ashres)
    post_sd[,i] = get_psd(ashres)
    lfsr[,i] = get_lfsr(ashres)
    loglik = loglik + get_loglik(ashres) #return the sum of loglikelihoods
  }
  posterior_matrices = list(PosteriorMean = post_mean,
                            PosteriorSD = post_sd,
                            lfsr = lfsr)

  m = list(result=posterior_matrices,loglik=loglik)
  class(m) = "mash_1by1"
  return(m)
}



#' Compute loglikelihood from a matrix of log-likelihoods and fitted pi
#' @param pi the vector of mixture proportions
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' @export
compute_loglik_from_matrix_and_pi = function(pi,lm){
  return(sum(log(lm$lik_matrix %*% pi)+lm$lfactors))
}


#' Initialize mixture proportions - currently by making them all equal
#' @param K the number of components
#' @return a vector of length K whose elements are positive and sums to 1
initialize_pi = function(K){
  return(rep(1/K,K))
}

grid_min = function(Bhat,Shat){
  min(Shat)/10
}

grid_max = function(Bhat,Shat){
  if (all(Bhat^2 <= Shat^2)) {
    8 * grid_min(Bhat,Shat) # the unusual case where we don't need much grid
  }  else {
    2 * sqrt(max(Bhat^2 - Shat^2))
  }
}

#' Automatically select grid
autoselect_grid = function(data,mult){
  gmax = grid_max(data$Bhat, data$Shat)
  gmin = grid_min(data$Bhat, data$Shat)
  if (mult == 0) {
    return(c(0, gmax/2))
  }
  else {
    npoint = ceiling(log2(gmax/gmin)/log2(mult))
    return(mult^((-npoint):0) * gmax)
  }
  #message("autoselect_grid is a place-holder\n")
  #return(c(0.5,1,2))
}

