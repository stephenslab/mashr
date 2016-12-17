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
#' @return post.means JxPxR array of posterior means, correspodning to the posterior mean ##' for the Jth individual in the Kth compoenent across all R tissues
#' @return post.covs JxPxR array of posterior vars, correspodning to the posterior vars ##' for the Jth individual in the Kth compoenent across all R tissues
#' @return post.ups JxPxR array of posterior tailprobs, corresponding to the marginal
#' @return upper tail probability for the Jth individual in the Pth compoenent across all R
#' @return post.downs JxKxR array of posterior tailprobs, corresponding to the marginal
#' @return lower  tail probability for the Jth individual in the Pth component across all R
#' @return post.nulls JxKxR array of posterior nullprobs, corresponding to the marginal
#' @return "null probability"" for the Jth individual in the Pth component across all R
#' @export

compute_posterior_arrays=function(Bhat,Shat,Ulist){
  R=ncol(Bhat)
  J=nrow(Bhat)
  P=length(Ulist)
  post_means=array(NA,dim=c(J,P,R))
  post_vars=array(NA,dim=c(J,P,R))
  post_ups=array(NA,dim=c(J,P,R))
  post_nulls=array(NA,dim=c(J,P,R))
  post_downs=array(NA,dim=c(J,P,R))

  for(j in 1:J){
    bhat=as.vector(t(Bhat[j,]))##turn i into a R x 1 vector
    Vhat=diag(Shat[j,])^2
    Vinv <- solve(Vhat)
    for(p in 1:P){
      U1 <- posterior_cov(Vinv, Ulist[[p]])
      mu1 <- as.array(posterior_mean(bhat, Vinv, U1))
      post_means[j,p,]=mu1
      post_vars[j,p,]=as.array(diag(U1))
      for(r in 1:R){##` marginal calculation for each tissue in each component
        if(post_vars[j,p,r]==0){###if the covariance matrix has entry 0, then p(null=1)
          post_ups[j,p,r]=0#need to figure out how to make the list have individual entries
          post_downs[j,p,r]=0
          post_nulls[j,p,r]=1}
        else{
          post_ups[j,p,r]=pnorm(0,mean=mu1[r],sd=sqrt(diag(U1)[r]),lower.tail=F)
          post_downs[j,p,r]=pnorm(0,mean=mu1[r],sd=sqrt(diag(U1)[r]),lower.tail=T)
          post_nulls[j,p,r]=0}
      }
    }
  }
  return(list(post_means=post_means,post_nulls=post_nulls,
              post_vars=post_vars,post_ups=post_ups,post_downs=post_downs))
}

