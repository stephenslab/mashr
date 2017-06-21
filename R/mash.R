#' Apply mash method to data
#' @param data a mash data object containing the Bhat matrix and standard errors; created using \code{set_mash_data}
#' @param Ulist a list of covariance matrices to use
#' @param gridmult scalar indicating factor by which adjacent grid values should differ; close to 1 for fine grid
#' @param grid vector of grid values to use (scaling factors omega in paper)
#' @param normalizeU whether or not to normalize the U covariances to have maximum of 1 on diagonal
#' @param usepointmass whether to include a point mass at 0, corresponding to null in every condition
#' @param g the value of g obtained from a previous mash fit - an alternative to supplying Ulist, grid and usepointmass
#' @param fixg if g is supplied, allows the mixture proportions to be fixed rather than estimated - e.g. useful for fitting mash to test data after fitting it to training data
#' @param prior indicates what penalty to use on the likelihood, if any
#' @param optmethod name of optimization method to use
#' @param verbose If \code{TRUE}, print progress to R console.
#' @param add.mem.profile If \code{TRUE}, print memory usage to R console (requires R library `profmem`).
#' @param algorithm.version Indicates whether to use R or Rcpp version
#' @param pi_thresh threshold below which mixture components are ignored in computing posterior summaries (to speed calculations by ignoring negligible components)
#' @return a list with elements result, loglik and fitted_g
#' @examples
#' Bhat = matrix(rnorm(100),ncol=5) # create some simulated data
#' Shat = matrix(rep(1,100),ncol=5)
#' data = mashr::set_mash_data(Bhat,Shat)
#' U.c = mashr::cov_canonical(data)
#' res.mash = mashr::mash(data,U.c)
#'
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
                optmethod = c("mixIP","mixEM","cxxMixSquarem"),
                verbose = TRUE,
                add.mem.profile = FALSE,
                algorithm.version = c("Rcpp","R"),
                pi_thresh = 1e-10) {

  algorithm.version = match.arg(algorithm.version)

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
  tryCatch(chol(data$V), error = function(e) stop("Input matrix V is not positive definite"))
  xUlist = expand_cov(Ulist,grid,usepointmass)

  # Check "add.mem.profile" argument.
  if (add.mem.profile)
    if (!requireNamespace("profmem",quietly = TRUE))
      stop("add.mem.profile = TRUE requires the profmem package")

  # Get the number of samples (J), the number of mixture components
  # (i.e., prior covariances).
  J <- nrow(data$Bhat)
  P <- length(xUlist)

  # Calculate likelihood matrix.
  if (verbose)
    cat(sprintf(" - Computing %d x %d likelihood matrix.\n",J,P))
  if (add.mem.profile) {
    out.time <- system.time(out.mem <- profmem::profmem({
      # lm <- calc_lik_matrix(data,xUlist,algorithm.version)
      lm <- calc_relative_lik_matrix(data,xUlist,log = TRUE,algorithm.version)
    },threshold = 1000))
  } else {
    out.time <-
      system.time(
    # lm<-calc_lik_matrix(data,xUlist,log=TRUE,algorithm.version)
    lm <- calc_relative_lik_matrix(data,xUlist,algorithm.version))
  }
  if (verbose) {
    if (add.mem.profile)
      cat(sprintf(paste(" - Likelihood calculations allocated %0.2f MB",
                        "and took %0.2f seconds.\n"),
                  sum(out.mem$bytes,na.rm = TRUE)/1024^2,
                  out.time["elapsed"]))
    else
      cat(sprintf(" - Likelihood calculations took %0.2f seconds.\n",
                  out.time["elapsed"]))
  }
    
  # Main fitting procedure.
  if(!fixg){
    if (verbose)
      cat(sprintf(" - Fitting model with %d mixture components.\n",P))
    prior <- set_prior(ncol(lm$lik_matrix),prior)
    if (add.mem.profile)
      out.time <- system.time(out.mem <- profmem::profmem({
        pi_s <- optimize_pi(lm$lik_matrix,prior=prior,optmethod=optmethod)
      },threshold = 1000))
    else
      out.time <- system.time(pi_s <-
                    optimize_pi(lm$lik_matrix,prior=prior,optmethod=optmethod))
    if (verbose)
      if (add.mem.profile)
        cat(sprintf(" - Model fitting allocated %0.2f MB and took %0.2f s.\n",
                     sum(out.mem$bytes,na.rm = TRUE)/1024^2,
                    out.time["elapsed"]))
      else
        cat(sprintf(" - Model fitting took %0.2f seconds.\n",
                    out.time["elapsed"]))
  }
  else{ #if fixg, just use g$pi for pi
    pi_s = g$pi
  }

  # threshold mixture components
  which.comp = (pi_s > pi_thresh)

  # Compute posterior matrices.
  posterior_weights <- compute_posterior_weights(pi_s[which.comp],lm$lik_matrix[,which.comp])
  if (verbose)
    cat(" - Computing posterior matrices.\n")
  if (add.mem.profile)
    out.time <- system.time(out.mem <- profmem::profmem({
      posterior_matrices <- compute_posterior_matrices(data,xUlist[which.comp],
                                                       posterior_weights, algorithm.version)
    },threshold = 1000))
  else
    out.time <-
      system.time(posterior_matrices <-
        compute_posterior_matrices(data,xUlist[which.comp],posterior_weights, algorithm.version))
  if (verbose)
    if (add.mem.profile)
      cat(sprintf(" - Computation allocated %0.2f MB and took %0.2f s.\n",
                  sum(out.mem$bytes,na.rm = TRUE)/1024^2,
                  out.time["elapsed"]))
    else
      cat(sprintf(" - Computation allocated took %0.2f seconds.\n",
                  out.time["elapsed"]))

  # Compute marginal log-likelihood.
  loglik = compute_loglik_from_matrix_and_pi(pi_s,lm)
  fitted_g = list(pi = pi_s, Ulist=Ulist, grid=grid, usepointmass=usepointmass)

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


# Sets prior to be a vector of length K depending on character string
# prior can be "nullbiased" or "uniform".
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
#'
#' @param data A list with the following two elements: \code{Bhat} an
#' n by R matrix of observations (n units in R conditions); and
#' \code{Shat}, an n by R matrix of standard errors (n units in R
#' conditions),
#'
#' @description Performs simple "condition-by-condition" analysis by
#' running \code{ash} from package \code{ashr} on data from each
#' condition, one at a time. May be a useful first step to identify
#' top hits in each condition before a mash analysis.
#'
#' @return A list similar to the output of mash, particularly
#' including posterior matrices.
#'
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
#' @param pi_s the vector of mixture proportions
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' @export
compute_loglik_from_matrix_and_pi = function(pi_s,lm){
  return(sum(log(lm$lik_matrix %*% pi_s)+lm$lfactors))
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

# Automatically select grid.
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

