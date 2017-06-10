#' @title calc_lik_vector
#' @description computes vector of likelihoods for bhat for each of P prior covariances
#' @param bhat Rx1 vector of bhat
#' @param V RxR covariance matrix for likelihood
#' @param Ulist list of prior covariance matrices
#' @param log if true computes log-likelihood
#' @return P vector of multivariate normal likelihoods, with pth element p(bhat | Ulist[p], V)
#' @importFrom mvtnorm dmvnorm
#' @export
calc_lik_vector=function(bhat,V,Ulist,log=FALSE){
  sapply(seq(1:length(Ulist)),function(p){dmvnorm(x=bhat, sigma=Ulist[[p]] + V,log=log)})
}

#' @title calc_lik_matrix
#' @description computes matrix of likelihoods for each of J rows of Bhat for each of P prior covariances
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param Ulist list of prior covariance matrices
#' @param log if true computes log-likelihood
#' @return J x P vector of multivariate normal likelihoods, p(bhat | Ulist[p], V)
#' @export
calc_lik_matrix = function(data, Ulist, log=FALSE){

  # TO DO: Implement a faster vresion of this function using Rcpp or
  # RcppArmadillo.
  J = n_effects(data)
  t(sapply(seq(1:J),
           function(j){
             calc_lik_vector(data$Bhat[j,],get_cov(data,j),Ulist,log)
             }))
}

#' @title Calculate matrix of relative likelihoods (likelihoods, normalized to avoid numeric issues)
#' @description computes matrix of relative likelihoods for each of J rows of Bhat for each of P prior covariances
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param Ulist list of prior covariance matrices
#' @return J x P matrix of likelihoods, p(bhat[j] | Ulist[p], V), but normalized so that the max in each row is 1
#' @export
calc_relative_lik_matrix = function(data, Ulist){
  matrix_llik = calc_lik_matrix(data,Ulist,log=TRUE)
  lfactors = apply(matrix_llik,1, max)
  matrix_llik = matrix_llik - lfactors #avoid numerical issues by subtracting max of each row
  return(list(lik_matrix = exp(matrix_llik), lfactors = lfactors ))
}

#' @title Calculate matrix of relative likelihoods C++ version 
#' @description computes matrix of relative likelihoods for each of J rows of Bhat for each of P prior covariances
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param Ulist list of prior covariance matrices
#' @return J x P matrix of likelihoods, p(bhat[j] | Ulist[p], V), but normalized so that the max in each row is 1
#' @useDynLib mashr
#' @exportPattern ^[[:alpha:]]+
#' @importFrom Rcpp evalCpp
#' @export
calc_relative_lik_matrix_arma = function(data, Ulist){
  matrix_llik = calc_lik_rcpp(t(data$Bhat), t(data$Shat), data$V, simplify2array(Ulist),log=TRUE)$data
  lfactors = apply(matrix_llik,1, max)
  matrix_llik = matrix_llik - lfactors #avoid numerical issues by subtracting max of each row
  return(list(lik_matrix = exp(matrix_llik), lfactors = lfactors ))
}
