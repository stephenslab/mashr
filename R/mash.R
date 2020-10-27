#' @title Apply mash method to data
#'
#' @param data a mash data object containing the Bhat matrix, standard
#' errors, alpha value; created using \code{mash_set_data} or
#' \code{mash_set_data_contrast}
#'
#' @param Ulist a list of covariance matrices to use 
#' (see \code{normalizeU} for rescaling these matrices) 
#'
#' @param gridmult scalar indicating factor by which adjacent grid
#' values should differ; close to 1 for fine grid
#'
#' @param grid vector of grid values to use (scaling factors omega in
#' paper)
#'
#' @param normalizeU whether or not to normalize the U covariances to
#' have maximum of 1 on diagonal
#'
#' @param usepointmass whether to include a point mass at 0,
#' corresponding to null in every condition
#'
#' @param g the value of g obtained from a previous mash fit - an
#' alternative to supplying Ulist, grid and usepointmass
#'
#' @param fixg if g is supplied, allows the mixture proportions to be
#' fixed rather than estimated; e.g., useful for fitting mash to test
#' data after fitting it to training data
#'
#' @param prior indicates what penalty to use on the likelihood, if any
#'
#' @param optmethod name of optimization method to use
#'
#' @param control A list of control parameters passed to optmethod.
#'
#' @param verbose If \code{TRUE}, print progress to R console.
#'
#' @param add.mem.profile If \code{TRUE}, print memory usage to R
#' console (requires R library `profmem`).
#'
#' @param algorithm.version Indicates whether to use R or Rcpp version
#'
#' @param pi_thresh threshold below which mixture components are
#' ignored in computing posterior summaries (to speed calculations by
#' ignoring negligible components)
#'
#' @param A the linear transformation matrix, Q x R matrix. This is
#' used to compute the posterior for Ab.
#'
#' @param posterior_samples the number of samples to be drawn from the
#' posterior distribution of each effect.
#'
#' @param seed A random number seed to use when sampling from the
#' posteriors. It is used when \code{posterior_samples > 0}.
#'
#' @param outputlevel controls amount of computation / output; 1:
#' output only estimated mixture component proportions, 2: and
#' posterior estimates, 3: and posterior covariance matrices, 4: and
#' likelihood matrices
#'
#' @return a list with elements result, loglik and fitted_g
#'
#' @examples
#' Bhat     = matrix(rnorm(100),ncol=5) # create some simulated data
#' Shat     = matrix(rep(1,100),ncol=5)
#' data     = mash_set_data(Bhat,Shat, alpha=1)
#' U.c      = cov_canonical(data)
#' res.mash = mash(data,U.c)
#'
#' @export
#'
mash = function(data,
                Ulist = NULL,
                gridmult= sqrt(2),
                grid = NULL,
                normalizeU = TRUE,
                usepointmass = TRUE,
                g = NULL,
                fixg = FALSE,
                prior=c("nullbiased","uniform"),
                optmethod = c("mixSQP","mixIP","mixEM","cxxMixSquarem"),
                control = list(),
                verbose = TRUE,
                add.mem.profile = FALSE,
                algorithm.version = c("Rcpp","R"),
                pi_thresh = 1e-10,
                A = NULL,
                posterior_samples = 0,
                seed = 123,
                outputlevel = 2) {

  algorithm.version = match.arg(algorithm.version)

  if(!missing(g)){ # g is supplied
    if(!missing(Ulist)){stop("cannot supply both g and Ulist")}
    if(!missing(grid)){stop("cannot supply both g and grid")}
    if(!missing(normalizeU)){stop("cannot supply both g and normalizeU")}
    if(!missing(usepointmass)){stop("cannot supply both g and usepointmass")}
    Ulist = g$Ulist
    grid = g$grid
    usepointmass = g$usepointmass
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

  # Get the number of samples (J), the number of mixture components
  # (i.e., prior covariances).
  J <- nrow(data$Bhat)
  R <- ncol(data$Bhat)

  for (i in 1:length(Ulist)) {
    check_covmat_basics(Ulist[[i]])
    if (nrow(Ulist[[i]]) != R)
      stop(paste("Matrices in Ulist must be of dimension", R, "by", R))
  }

  xUlist = expand_cov(Ulist,grid,usepointmass)
  P <- length(xUlist)

  # Check "add.mem.profile" argument.
  if (add.mem.profile)
    if (!requireNamespace("profmem",quietly = TRUE))
      stop("add.mem.profile = TRUE requires the profmem package")


  # Calculate likelihood matrix.
  if (verbose)
    cat(sprintf(" - Computing %d x %d likelihood matrix.\n",J,P))
  if (add.mem.profile) {
    out.time <- system.time(out.mem <- profmem::profmem({
      lm <- calc_relative_lik_matrix(data,xUlist, algorithm.version)
    },threshold = 1000))
  } else {
    out.time <- system.time(
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
    prior <- set_prior(ncol(lm$loglik_matrix),prior)
    if (add.mem.profile)
      out.time <- system.time(out.mem <- profmem::profmem({
        pi_s <- optimize_pi(exp(lm$loglik_matrix),prior=prior,optmethod=optmethod, control=control)
      },threshold = 1000))
    else
      out.time <- system.time(pi_s <-
                    optimize_pi(exp(lm$loglik_matrix),prior=prior,optmethod=optmethod, control=control))
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
  posterior_weights <- compute_posterior_weights(pi_s[which.comp],exp(lm$loglik_matrix[,which.comp]))
  # Compute posterior matrices.
  if (outputlevel > 1) {
    if (verbose)
      cat(" - Computing posterior matrices.\n")
    if (add.mem.profile)
      out.time <- system.time(out.mem <- profmem::profmem({
        posterior_matrices <- compute_posterior_matrices(data, xUlist[which.comp],
                                                         posterior_weights, algorithm.version, A=A,
                                                         output_posterior_cov=(outputlevel > 2),
                                                         posterior_samples = posterior_samples, seed = seed)
      },threshold = 1000))
    else
      out.time <-
        system.time(posterior_matrices <-
          compute_posterior_matrices(data,xUlist[which.comp],
                                     posterior_weights, algorithm.version, A=A,
                                     output_posterior_cov=(outputlevel > 2),
                                     posterior_samples = posterior_samples, seed = seed))
    if (verbose)
      if (add.mem.profile)
        cat(sprintf(" - Computation allocated %0.2f MB and took %0.2f s.\n",
                    sum(out.mem$bytes,na.rm = TRUE)/1024^2,
                    out.time["elapsed"]))
      else
        cat(sprintf(" - Computation allocated took %0.2f seconds.\n",
                    out.time["elapsed"]))
  } else {
    posterior_matrices = NULL
  }
  # Compute marginal log-likelihood.
  vloglik = compute_vloglik_from_matrix_and_pi(pi_s,lm,data$Shat_alpha)
  loglik = sum(vloglik)
  if(usepointmass){ # compute BF
    null_loglik = compute_null_loglik_from_matrix(lm,data$Shat_alpha)
    alt_loglik = compute_alt_loglik_from_matrix_and_pi(pi_s,lm,data$Shat_alpha)
  } else {
    null_loglik = NULL
    alt_loglik = NULL
  }
  # results
  fitted_g = list(pi=pi_s, Ulist=Ulist, grid=grid, usepointmass=usepointmass)
  names(posterior_weights) = which(which.comp)
  m=list(result = posterior_matrices,
         loglik = loglik, vloglik = vloglik,
         null_loglik = null_loglik,
         alt_loglik = alt_loglik,
         fitted_g = fitted_g,
         posterior_weights = posterior_weights,
         alpha=data$alpha)
  #for debugging
  if(outputlevel==4){m = c(m,list(lm=lm))}
  class(m) = "mash"
  return(m)
}

#' @title Compute posterior matrices for fitted mash object on new
#' data
#'
#' @param g a mash object or the fitted_g from a mash object.
#'
#' @param data a set of data on which to compute the posterior
#' matrices
#'
#' @param pi_thresh threshold below which mixture components are
#' ignored in computing posterior summaries (to speed calculations by
#' ignoring negligible components)
#'
#' @param algorithm.version Indicates whether to use R or Rcpp version
#'
#' @param A the linear transformation matrix, Q x R matrix. This is
#' used to compute the posterior for Ab.
#'
#' @param posterior_samples the number of samples to be drawn from the
#' posterior distribution of each effect.
#'
#' @param seed a random number seed to use when sampling from the
#' posteriors. It is used when \code{posterior_samples > 0}.
#'
#' @param output_posterior_cov whether or not to output posterior
#' covariance matrices for all effects
#'
#' @return A list of posterior matrices
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data))
#' mash_compute_posterior_matrices(m,data)
#'
#' @export
#'
mash_compute_posterior_matrices = function(g, data, pi_thresh = 1e-10, algorithm.version = c("Rcpp", "R"), A=NULL, output_posterior_cov=FALSE,
                                           posterior_samples = 0, seed = 123){

  if(class(g)=="mash"){
    alpha = g$alpha
    g = g$fitted_g
    if(alpha != data$alpha){
      stop('The alpha in data is not the one used to compute the mash model.')
    }
  }
  else{
    message('Warning: Please make sure the alpha in data is consistent with the `alpha` used to compute the fitted_g.')
  }

  xUlist = expand_cov(g$Ulist,g$grid,g$usepointmass)
  lm_res = calc_relative_lik_matrix(data, xUlist, algorithm.version = algorithm.version)
  which.comp = (g$pi > pi_thresh)

  print(system.time(
  posterior_weights <- compute_posterior_weights(g$pi[which.comp], exp(lm_res$loglik_matrix[,which.comp]))))
  print(system.time(
  posterior_matrices <- compute_posterior_matrices(data, xUlist[which.comp],
                                                  posterior_weights,
                                                  algorithm.version, A=A, output_posterior_cov=output_posterior_cov,
                                                  posterior_samples = posterior_samples, seed=seed)))
  names(posterior_weights) = which(which.comp)
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

#' @title Create expanded list of covariance matrices expanded by
#'   grid, Sigma_{lk} = omega_l U_k
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param Ulist a list of covarance matrices
#'
#' @param grid a grid of scalar values by which the covariance
#'   matrices are to be sc
#'
#' @param usepointmass if TRUE adds a point mass at 0 (null component)
#'   to the list
#'
#' @return This takes the covariance matrices in Ulist and multiplies
#' them by the grid values If usepointmass is TRUE then it adds a null
#' component.
#'
#' @keywords internal
#'
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
#' @param alpha Numeric value of alpha parameter in the model. alpha =
#' 0 for Exchangeable Effects (EE), alpha = 1 for Exchangeable
#' Z-scores (EZ).
#'
#' @param ... optionally, other parameters to be passed to ash
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
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' mash_1by1(simdata)
#'
#' @export
#'
mash_1by1 = function(data, alpha=0, ...){
  Bhat = data$Bhat
  Shat = data$Shat
  post_mean = post_sd = lfsr = matrix(nrow = nrow(Bhat), ncol= ncol(Bhat))
  loglik = 0
  for(i in 1:ncol(Bhat)){
    ashres = ash(Bhat[,i],Shat[,i],mixcompdist="normal", alpha=alpha, ...)
    post_mean[,i] = get_pm(ashres)
    post_sd[,i] = get_psd(ashres)
    lfsr[,i] = get_lfsr(ashres)
    loglik = loglik + get_loglik(ashres) #return the sum of loglikelihoods across conditions
  }
  posterior_matrices = list(PosteriorMean = post_mean,
                            PosteriorSD = post_sd,
                            lfsr = lfsr)
  for (i in names(posterior_matrices)) {
    if (!is.null(colnames(data$Bhat))) colnames(posterior_matrices[[i]]) = colnames(data$Bhat)
    if (!is.null(rownames(data$Bhat))) rownames(posterior_matrices[[i]]) = rownames(data$Bhat)
  }

  m = list(result=posterior_matrices,loglik=loglik)
  class(m) = "mash_1by1"
  return(m)
}

# @title Initialize mixture proportions - currently by making them all equal
# @param K the number of components
# @return a vector of length K whose elements are positive and sums to 1
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
  include = !(data$Shat==0 | !is.finite(data$Shat) | is.na(data$Bhat))
  gmax = grid_max(data$Bhat[include], data$Shat[include])
  gmin = grid_min(data$Bhat[include], data$Shat[include])
  if (mult == 0) {
    return(c(0, gmax/2))
  }
  else {
    npoint = ceiling(log2(gmax/gmin)/log2(mult))
    return(mult^((-npoint):0) * gmax)
  }
}
