#' @title Compute vector of loglikelihood for fitted mash object on
#'   new data
#'
#' @param g A mash object.
#'
#' @param data A set of data on which to compute the loglikelihood.
#'
#' @param algorithm.version Indicate R or Rcpp version
#'
#' @return The vector of log-likelihoods for each data point computed
#' using g.
#'
#' @details The log-likelihood for each element is \eqn{p(Bhat_j |
#' Shat_j,g,\alpha)} where \eqn{Bhat_j | B_j, Shat_j \sim N(B_j,
#' Shat_j)} and \eqn{B_j/Shat_j^\alpha | Shat_j \sim g} Here the value
#' of \eqn{\alpha} is set when setting up the data object in
#' `mash_set_data`. If g is a mash object (safest!) then the function
#' will check that this value matches the \eqn{\alpha} used when
#' fitting `mash`. Note: as a convenience, this function can also be
#' called with g a mixture distribution with same structure as the
#' fitted_g from a mash object. This is mostly useful when doing
#' simulations, where you might want to compute the likelihood under
#' the "true" g. When used in this way the user is responsible for
#' making sure that the g makes sense with the alpha set in data.
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data))
#' mash_compute_vloglik(m,data)
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
#'   Shat_j,g,\alpha)} where \eqn{Bhat_j | B_j, Shat_j \sim N(B_j,
#'   Shat_j)} and \eqn{B_j/Shat_j^\alpha | Shat_j \sim g}.
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data))
#' mash_compute_loglik(m,data)
#'
#' @export
#'
mash_compute_loglik = function(g,data, algorithm.version=c("Rcpp", "R")){
  return( sum( mash_compute_vloglik(g,data, algorithm.version = algorithm.version) ) )
}

#' @title Compute vector of alternative loglikelihoods from a matrix
#'   of log-likelihoods and fitted pi
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param pi_s the vector of mixture proportions, with first element
#' corresponding to null
#'
#' @param lm the results of a likelihood matrix calculation from
#'   \code{calc_relative_lik_matrix} whose first column corresponds to
#'   null
#'
#' @param Shat_alpha matrix of Shat^alpha
#'
#' @keywords internal
#'
compute_alt_loglik_from_matrix_and_pi = function(pi_s,lm,Shat_alpha){
  if(pi_s[1] == 1){
    # 1-pi_s[1] = 0 --> NaN
    # so we use weight 1/(P-1), where P is the length of pi_s
    tmp = rep(1/(length(pi_s)-1), (length(pi_s)-1))
    return(log(exp(lm$loglik_matrix[,-1,drop=FALSE]) %*% (tmp))+lm$lfactors-rowSums(log(Shat_alpha)))
  }else{
    return(log(exp(lm$loglik_matrix[,-1,drop=FALSE]) %*% (pi_s[-1]/(1-pi_s[1])))+lm$lfactors-rowSums(log(Shat_alpha)))
  }
}

#' @title Compute a vector of null loglikelihoods from a matrix of
#'   log-likelihoods
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param lm the results of a likelihood matrix calculation from
#'   \code{calc_relative_lik_matrix} whose first column corresponds to
#'   null
#'
#' @param Shat_alpha matrix of Shat^alpha
#'
#' @keywords internal
#'
compute_null_loglik_from_matrix = function(lm,Shat_alpha){
  lm$loglik_matrix[,1]+lm$lfactors-rowSums(log(Shat_alpha))
}

#' @title Computes a vector of loglikelihoods from a matrix of
#'   log-likelihoods and fitted pi
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param pi_s the vector of mixture proportions
#'
#' @param lm the results of a likelihood matrix calculation from
#'   \code{calc_relative_lik_matrix}
#'
#' @param Shat_alpha matrix of Shat^alpha
#'
#' @keywords internal
#'
compute_vloglik_from_matrix_and_pi = function(pi_s,lm,Shat_alpha){
  return(log(exp(lm$loglik_matrix) %*% pi_s)+lm$lfactors-rowSums(log(Shat_alpha)))
}

#' @title Compute the total loglikelihood from a matrix of
#'   log-likelihoods and fitted pi
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param pi_s the vector of mixture proportions
#'
#' @param lm the results of a likelihood matrix calculation from
#'   \code{calc_relative_lik_matrix}
#'
#' @param Shat_alpha matrix of Shat^alpha
#'
#' @keywords internal
#'
compute_loglik_from_matrix_and_pi = function(pi_s,lm,Shat_alpha){
  return(sum(compute_vloglik_from_matrix_and_pi(pi_s,lm,Shat_alpha)))
}
