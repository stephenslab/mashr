#' @title calc_lik_matrix_common_cov
#' @description computes matrix of likelihoods for each of J rows of Bhat for each of P prior covariances; special case when standard errors and variances are all same across j;
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param Ulist list of prior covariance matrices
#' @param log if true computes log-likelihood
#' @return J x P vector of multivariate normal likelihoods, p(bhat | Ulist[p], V), where V is same for each bhat
#' @details Compared with \code{calc_lik_matrix} this function exploits fact that observations are iid in this case,
#' so the inverse covariance matrices only need to be done once, reducing computation to R^3 + JR^2 instead of JR^3
#' @export
calc_lik_matrix_common_cov = function(data, Ulist, log=FALSE){
  V = mashr:::get_cov(data,1) # all covariances are same
  res = plyr::laply(Ulist,
    function(U){dmvnorm(x=data$Bhat,sigma=V+U,log=log)})
  dimnames(res) = NULL # just to make result identical to the non-common-cov version
  t(res)
}

#' @title Check that all covariates are equal.
#'
#' @description checks if all rows of Shat are the same - if so
#'     covariances are equal
#'
#' @param data A mash data object.
#'
#' @export
is_common_cov = function(data){
  all((t(data$Shat) - data$Shat[1,]) == 0)
}
