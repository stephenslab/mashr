#' @title calc_lik_vector
#' @details computes vector of likelihoods for bhat for each of P prior covariances
#' @param bhat Rx1 vector of bhat
#' @param V RxR covariance matrix for likelihood
#' @param Ulist list of prior covariance matrices
#' @param log if true computes log-likelihood
#' @return P vector of multivariate normal likelihoods, with pth element p(bhat | Ulist[p], V)
#' @export
calc_lik_vector=function(bhat,V,Ulist,log=FALSE){
  sapply(seq(1:length(Ulist)),function(p){mvtnorm::dmvnorm(x=bhat, sigma=Ulist[[p]] + V,log=log)})
}

#' @title calc_lik_matrix
#' @details computes matrix of likelihoods for each of J rows of Bhat for each of P prior covariances
#' @param Bhat J x R matrix of betahat values
#' @param Shat J x R matrix of standard errors
#' @param Ulist list of prior covariance matrices
#' @param log if true computes log-likelihood
#' @return J x P vector of multivariate normal likelihoods, p(bhat | Ulist[p], V)
#' @export
calc_lik_matrix = function(Bhat, Shat, Ulist, log=FALSE){
  J = nrow(Bhat)
  t(sapply(seq(1:J),function(j){calc_lik_vector(Bhat[j,],diag(Shat[j,]^2),Ulist,log)}))
}

