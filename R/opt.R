#' @title Estimate mixture weights by maximum (penalized) likelihood
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param matrix_lik a matrix of likelihoods, where the (i,k)th entry
#'   is the probability of observation i given it came from component k
#'   of g
#'
#' @param pi_init numeric vector specifying value from which to
#'   initialize optimization
#'
#' @param prior numeric vector specifying prior to use in the
#'   penalized likelihood
#'
#' @param optmethod a string, giving name of optimization function to
#'   use
#'
#' @param control a list of parameters to be passed to optmethod
#'
#' @return numeric vector specifying the optimal mixture weights
#'
#' @importFrom assertthat are_equal
#' @importFrom ashr mixIP
#' @importFrom ashr mixEM
#' @importFrom ashr cxxMixSquarem
#' @importFrom ashr mixSQP
#'
#' @keywords internal
#'
optimize_pi = function(matrix_lik, pi_init = NULL,
                       prior=NULL,
                       optmethod=c("mixSQP","mixIP","mixEM","cxxMixSquarem"),
                       control=list() ){
  optmethod = match.arg(optmethod)

  if (optmethod == "mixIP")
   if (!requireNamespace("REBayes",quietly = TRUE)) {
     warning(paste("optmethod = \"mixIP\" requires REBayes package;",
                   "switching to optmethod = \"mixSQP\""))
     optmethod <- "mixSQP"
   }

  K = ncol(matrix_lik)
  if(missing(prior)){prior = rep(1,K)}
  if(missing(pi_init)){pi_init = initialize_pi(K)}
  are_equal(length(pi_init),K)
  are_equal(length(prior),K)
  res = do.call(optmethod, args= list(matrix_lik = matrix_lik, prior=prior, pi_init = pi_init,control=control))
  return(res$pihat)
}

