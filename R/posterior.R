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

#' @title compute_posterior_arrays
#' @param Bhat J by R matrix of estimates
#' @param Shat J x R matrix of standard errors
#' @param Ulist list of P prior covariance matrices
#' @return post_means JxPxR array of posterior means
#' @return post_vars JxPxR array of posterior (marginal) variances
#' @return post_pos JxPxR array of posterior (marginal) probability of being positive
#' @return post_neg JxPxR array of posterior (marginal) probability of being negative
#' @return post_zero JxPxR array of posterior (marginal) probability of being zero
#' @export

compute_posterior_arrays=function(Bhat,Shat,Ulist){
  R=ncol(Bhat)
  J=nrow(Bhat)
  P=length(Ulist)
  post_mean=array(NA,dim=c(J,P,R))
  post_var=array(NA,dim=c(J,P,R))
  post_pos=array(NA,dim=c(J,P,R))
  post_zero=array(NA,dim=c(J,P,R))
  post_neg=array(NA,dim=c(J,P,R))

  for(j in 1:J){
    bhat=as.vector(t(Bhat[j,]))##turn i into a R x 1 vector
    V=diag(Shat[j,])^2
    Vinv <- solve(V)
    for(p in 1:P){
      U1 <- posterior_cov(Vinv, Ulist[[p]])
      mu1 <- as.array(posterior_mean(bhat, Vinv, U1))
      post_mean[j,p,]= mu1
      post_var[j,p,]= diag(U1)
      post_pos[j,p,] = ifelse(post_var[j,p,]==0,0,pnorm(0,mean=mu1,sqrt(diag(U1)),lower.tail=F))
      post_neg[j,p,] = ifelse(post_var[j,p,]==0,0,pnorm(0,mean=mu1,sqrt(diag(U1)),lower.tail=T))
      post_zero[j,p,] = ifelse(post_var[j,p,]==0,1,0)
    }
  }
  return(list(post_mean=post_mean,post_zero=post_zero,
              post_var=post_var,post_pos=post_pos,post_neg=post_neg))
}

