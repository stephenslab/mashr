#' @title Compute vector of loglikelihood for fitted mash object on new data.
#'
#' @param g A mash object.
#'
#' @param data A set of data on which to compute the loglikelihood.
#'
#' @param algorithm.version Indicate R or Rcpp version
#'
#' @return The vector of log-likelihoods for each data point computed using g.
#'
#' @details The log-likelihood for each element is \eqn{p(Bhat_j |
#' Shat_j,g,\alpha)} where \eqn{Bhat_j | B_j, Shat_j \sim N(B_j,
#' Shat_j)} and \eqn{B_j/Shat_j^\alpha | Shat_j \sim g} Here the value
#' of \eqn{\alpha} is set when setting up the data object in
#' `mash_set_data`. If g is a mash object (safest!) then the function
#' will check that this value matches the \eqn{\alpha} used when
#' fitting `mash`.  Note: as a convenience, this function can also be
#' called with g a mixture distribution with same structure as the
#' fitted_g from a mash object. This is mostly useful when doing
#' simulations, where you might want to compute the likelihood under
#' the "true" g. When used in this way the user is responsible for
#' making sure that the g makes sense with the alpha set in data.
#'
#' @export
#'
mash_compute_vloglik = function(g,data, algorithm.version= c("Rcpp","R")){
  if(class(g)=="mash"){
    alpha = g$alpha
    g = g$fitted_g
    if(alpha != data$alpha){
      stop('The alpha in data does not match the one used to fit the mash model.')
    }
  }
  else{
    message('Warning: Please make sure the alpha in data is consistent with the `alpha` used to compute g.')
  }

  xUlist = expand_cov(g$Ulist,g$grid,g$usepointmass)
  lm_res = calc_relative_lik_matrix(data,xUlist,algorithm.version=algorithm.version)
  return(log(exp(lm_res$loglik_matrix) %*% g$pi) + lm_res$lfactors - rowSums(log(data$Shat_alpha)))
}

#' @title Compute loglikelihood for fitted mash object on new data.
#'
#' @param g A mash object or the fitted_g from a mash object.
#'
#' @param data A set of data on which to compute the loglikelihood.
#'
#' @param algorithm.version Indicate R or Rcpp version
#'
#' @return The log-likelihood for data computed using g.
#'
#' @details The log-likelihood for each element is \eqn{p(Bhat_j |
#' Shat_j,g,\alpha)} where \eqn{Bhat_j | B_j, Shat_j \sim N(B_j,
#' Shat_j)} and \eqn{B_j/Shat_j^\alpha | Shat_j \sim g}.
#'
#' @export
#'
mash_compute_loglik = function(g,data, algorithm.version=c("Rcpp", "R")){
  return( sum( mash_compute_vloglik(g,data, algorithm.version = algorithm.version) ) )
}

#' Compute vector of alternative loglikelihoods from a matrix of log-likelihoods and fitted pi
#' @param pi_s the vector of mixture proportions, with first element corresponding to null
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' whose first column corresponds to null
#' @param Shat_alpha matrix of Shat^alpha
compute_alt_loglik_from_matrix_and_pi = function(pi_s,lm,Shat_alpha){
  return(log(exp(lm$loglik_matrix[,-1,drop=FALSE]) %*% (pi_s[-1]/(1-pi_s[1])))+lm$lfactors-rowSums(log(Shat_alpha)))
}

#' Compute vector of null loglikelihoods from a matrix of log-likelihoods
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' whose first column corresponds to null
#' @param Shat_alpha matrix of Shat^alpha
compute_null_loglik_from_matrix = function(lm,Shat_alpha){
  lm$loglik_matrix[,1]+lm$lfactors-rowSums(log(Shat_alpha))
}

#' Compute vector of loglikelihoods from a matrix of log-likelihoods and fitted pi
#' @param pi_s the vector of mixture proportions
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' @param Shat_alpha matrix of Shat^alpha
compute_vloglik_from_matrix_and_pi = function(pi_s,lm,Shat_alpha){
  return(log(exp(lm$loglik_matrix) %*% pi_s)+lm$lfactors-rowSums(log(Shat_alpha)))
}

#' Compute total loglikelihood from a matrix of log-likelihoods and fitted pi
#' @param pi_s the vector of mixture proportions
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' @param Shat_alpha matrix of Shat^alpha
compute_loglik_from_matrix_and_pi = function(pi_s,lm,Shat_alpha){
  return(sum(compute_vloglik_from_matrix_and_pi(pi_s,lm,Shat_alpha)))
}
