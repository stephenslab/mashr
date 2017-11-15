#' Compute vector of loglikelihood for fitted mash object on new data
#' @param g a mash object
#' @param data a set of data on which to compute the loglikelihood
#' @return the vector of log-likelihoods for each data point computed using g
#' @details The log-likelihood for each element is $p(Bhat_j | Shat_j,g,\alpha)$
#' where $Bhat_j | B_j, Shat_j \sim N(B_j, Shat_j)$ and $B_j/Shat_j^\alpha | Shat_j \sim g$
#' Here the value of $\alpha$ is set when setting up the data object in `set_mash_data`.
#' If g is a mash object (safest!) then the function will check that this value matches the $\alpha$ used when fitting `mash`.
#' Note: as a convenience, this function can also be called with g a mixture distribution with same structure as the fitted_g from a mash object.
#' This is mostly useful when doing simulations, where you might want to compute the
#' likelihood under the "true" g. When used in this way the user is responsible for
#' making sure that the g makes sense with the alpha set in data.
#' @export
mash_compute_vloglik = function(g,data){
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
  lm_res = calc_relative_lik_matrix(data,xUlist)
  return(log(lm_res$lik_matrix %*% g$pi) + lm_res$lfactors - rowSums(log(data$Shat_alpha)))
}

#' Compute loglikelihood for fitted mash object on new data
#' @param g a mash object or the fitted_g from a mash object
#' @param data a set of data on which to compute the loglikelihood
#' @return the log-likelihood for data computed using g
#' @details The log-likelihood for each element is $p(Bhat_j | Shat_j,g,\alpha)$
#' where $Bhat_j | B_j, Shat_j \sim N(B_j, Shat_j)$ and $B_j/Shat_j^\alpha | Shat_j \sim g$.
#' @export
mash_compute_loglik = function(g,data){
  return( sum( mash_compute_vloglik(g,data) ) )
}

#' Compute vector of alternative loglikelihoods from a matrix of log-likelihoods and fitted pi
#' @param pi_s the vector of mixture proportions, with first element corresponding to null
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' whose first column corresponds to null
#' @param Shat_alpha matrix of Shat^alpha
compute_alt_loglik_from_matrix_and_pi = function(pi_s,lm,Shat_alpha){
  return(log(lm$lik_matrix[,-1,drop=FALSE] %*% (pi_s[-1]/(1-pi_s[1])))+lm$lfactors-rowSums(log(Shat_alpha)))
}

#' Compute vector of null loglikelihoods from a matrix of log-likelihoods
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' whose first column corresponds to null
#' @param Shat_alpha matrix of Shat^alpha
compute_null_loglik_from_matrix = function(lm,Shat_alpha){
  log(lm$lik_matrix[,1])+lm$lfactors-rowSums(log(Shat_alpha))
}

#' Compute vector of loglikelihoods from a matrix of log-likelihoods and fitted pi
#' @param pi_s the vector of mixture proportions
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' @param Shat_alpha matrix of Shat^alpha
compute_vloglik_from_matrix_and_pi = function(pi_s,lm,Shat_alpha){
  return(log(lm$lik_matrix %*% pi_s)+lm$lfactors-rowSums(log(Shat_alpha)))
}

#' Compute total loglikelihood from a matrix of log-likelihoods and fitted pi
#' @param pi_s the vector of mixture proportions
#' @param lm the results of a likelihood matrix calculation from \code{calc_relative_lik_matrix}
#' @param Shat_alpha matrix of Shat^alpha
compute_loglik_from_matrix_and_pi = function(pi_s,lm,Shat_alpha){
  return(sum(compute_vloglik_from_matrix_and_pi(pi_s,lm,Shat_alpha)))
}
