#' @title posterior_cov
#' @param Vinv R x R inverse covariance matrix for the likelihood
#' @param U R x R prior covariance matrix
#' @return R x R posterior covariance matrix
#' @details If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This function returns U1.
#' @export
posterior_cov <- function(Vinv, U){
  return(U %*% solve(Vinv %*% U + diag(nrow(U))))
}


#' @title posterior_mean
#' @param bhat R vector of observations
#' @param Vinv R x R inverse covariance matrix for the likelihood
#' @param U1 R x R posterior covariance matrix, computed using posterior_cov
#' @return R vector of posterior mean
#' @details If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This function returns mu1.
#' @export
posterior_mean <- function(bhat, Vinv, U1){
  return(U1 %*% Vinv %*% bhat)
}

#' @title Compute posterior matrices
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param Ulist a list of covariance matrices for each mixture component
#' @param posterior_weights the posterior probabilities of each mixture component in Ulist for the data
#' @return post_mean JxR matrix of posterior means
#' @return post_sd JxR matrix of posterior (marginal) standard deviations
#' @return post_pos JxR matrix of posterior (marginal) probability of being positive
#' @return post_neg JxR matrix of posterior (marginal) probability of being negative
#' @return post_zero JxR matrix of posterior (marginal) probability of being zero
#' @export
compute_posterior_matrices=function(data,Ulist,posterior_weights){

  post_arrays = compute_posterior_arrays(data,Ulist)
  post_mean=compute_weighted_quantity(post_arrays$post_mean,posterior_weights)
  post_mean2=compute_weighted_quantity(post_arrays$post_mean2,posterior_weights)
  post_sd = sqrt(post_mean2 - post_mean^2)
  #post_pos=compute_weighted_quantity(post_arrays$post_pos,posterior_weights)
  post_zero=compute_weighted_quantity(post_arrays$post_zero,posterior_weights)
  post_neg=compute_weighted_quantity(post_arrays$post_neg,posterior_weights)
  lfsr = ashr:::compute_lfsr(post_neg,post_zero)
  return(list(post_mean = post_mean,
              post_sd = post_sd,
              post_zero = post_zero,
              post_neg = post_neg,
              lfsr = lfsr))
}

#' @title compute_posterior_arrays
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param Ulist list of P prior covariance matrices
#' @return post_mean JxPxR array of posterior means
#' @return post_mean2 JxPxR array of posterior second moments
#' @return post_var JxPxR array of posterior variances
#' @return post_pos JxPxR array of posterior (marginal) probability of being positive
#' @return post_neg JxPxR array of posterior (marginal) probability of being negative
#' @return post_zero JxPxR array of posterior (marginal) probability of being zero
#' @export
compute_posterior_arrays=function(data,Ulist){
  R=n_conditions(data)
  J=n_effects(data)

  P=length(Ulist)
  post_mean=array(NA,dim=c(J,P,R))
  post_mean2 = array(NA,dim=c(J,P,R)) #mean squared value
  post_var = array(NA,dim=c(J,P,R)) #mean squared value
  post_pos=array(NA,dim=c(J,P,R))
  post_zero=array(NA,dim=c(J,P,R))
  post_neg=array(NA,dim=c(J,P,R))

  for(j in 1:J){
    bhat=as.vector(t(data$Bhat[j,]))##turn i into a R x 1 vector
    V=diag(data$Shat[j,])^2
    Vinv <- solve(V)
    for(p in 1:P){
      U1 <- posterior_cov(Vinv, Ulist[[p]])
      mu1 <- as.array(posterior_mean(bhat, Vinv, U1))
      post_mean[j,p,]= mu1
      post_var[j,p,] = diag(U1)
      post_mean2[j,p,] = mu1^2 + diag(U1) #diag(U1) is the posterior variance
      post_pos[j,p,] = ifelse(post_var[j,p,]==0,0,pnorm(0,mean=mu1,sqrt(diag(U1)),lower.tail=F))
      post_neg[j,p,] = ifelse(post_var[j,p,]==0,0,pnorm(0,mean=mu1,sqrt(diag(U1)),lower.tail=T))
      post_zero[j,p,] = ifelse(post_var[j,p,]==0,1,0)
    }
  }
  return(list(post_mean=post_mean,
              post_var = post_var,
              post_zero=post_zero,
              post_mean2= post_mean2,
              post_pos=post_pos,
              post_neg=post_neg))
}


#' @title Compute weighted means of posterior arrays
#' @description Generates a K x R matrix of posterior quantities (eg posterior mean) for each effect
#' @param post_array J x K x R array of posterior quantity for each effect for each component in each condition
#' @param weights J x K matrix of weights for each effect in each component (usually the posterior weights)
#' @return J by R matrix of quantities (eg posterior mean) for each effect in each condition. The (j,r) element is sum_k pi[j,k] a[j,k,r]
#' @export
compute_weighted_quantity = function(post_array,posterior_weights){
  R = dim(post_array)[3]
  weighted_array = post_array * rep(posterior_weights,R)
  return(apply(weighted_array,c(1,3),sum))
}

#' @title compute posterior probabilities
#' @details computes posterior probabilities that each effect came from each component
#' @param pi a K vector of mixture proportions
#' @param lik_mat a JxK matrix of likelihoods
#' @return a JxK matrix of posterior probabilities, the jth row contains posteriors for jth effect
compute_posterior_weights=function(pi,lik_mat){
  d= t(pi * t(lik_mat))
  norm = rowSums(d) # normalize probabilities to sum to 1
  return(d/norm)
}


