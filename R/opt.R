#' Estimate the mixture weights by maximum (penalized) likelihood
#' @param matrix_lik a matrix of likelihoods, where the (i,k)th entry is the probability of observation i given it came from component k of g
#' @param pi_init numeric vector specifying value from which to initialize optimization
#' @param prior numeric vector specifying prior to use in the penalized likelihood
#' @param optmethod a string, giving name of optimization function to use
#' @param control a list of parameters to be passed to optmethod
#' @return numeric vector specifying the optimal mixture weights
optimize_pi = function(matrix_lik, pi_init = NULL,
                       prior=NULL,
                       optmethod=c("mixIP","cxxMixSquarem"),
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

