#' @title Estimate null correlations
#' @description Estimates a null correlation matrix from data
#' @param data a mash data object, eg as created by \code{mash_set_data}
#' @param z_thresh the z score threshold below which to call an effect null
#' @param apply_lower_bound indicates whether to also apply a lower bound for the
#' correlations which can be computed based on all data. See \code{compute_null_correlation_lower_bound}
#' @details Returns the empirical correlation matrix of the effects that are "null" based on simple z score threshold
#' @importFrom stats cor
#' @export
estimate_null_correlation = function(data, z_thresh=2, apply_lower_bound = FALSE){
  z = data$Bhat/data$Shat
  max_absz = apply(abs(z),1, max)
  nullish = which(max_absz < z_thresh)
  if(length(nullish)<n_conditions(data)){
    stop("not enough null data to estimate null correlation")
  }
  nullish_z = z[nullish,]
  Vhat = cor(nullish_z)
  if(apply_lower_bound){
    lb = compute_null_correlation_lower_bound(data)
    Vhat = ifelse(Vhat>lb, Vhat, lb)
  }
  return(Vhat)
}

#' @title compute_null_correlation_lower_bound
#' @description Compute a lower bound on the null correlation matrix from data
#' @details Assume $$z^j_r = mu^j_r + e_r$$ where $e_1,e_2$ are joint normal with
#' variance 1 and some covariance (same as correlation since they are variance 1).
#' Then $E((z^j_1 - z^j_2)^2) = E(mu^j_1-mu^j_2)^2 + 2(1-cov(e_1,e_2)) > 2(1-cov(e_1,e_2))$.
#' Thus $$cov(e_1,e_2) > 1- 0.5E((z^j_1 - z^j_2)^2)$$ gives a lower bound on the covariance.
#'
#' @param data Description of this argument goes here.
compute_null_correlation_lower_bound = function(data){
  R = n_conditions(data)
  z = data$Bhat/data$Shat
  lb = matrix(0,ncol=R,nrow=R) # the lower bound
  for(i in 1:(R-1)){
    for(j in (i+1):R){
      lb[i,j] = 1-0.5*mean((z[,i]-z[,j])^2)
      lb[j,i] = lb[i,j] #lb is symmetric
    }
  }
  return(lb)
}
