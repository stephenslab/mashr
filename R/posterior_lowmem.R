#' @title Compute posterior matrices (low memory version)
#' @description Computes posterior matrices without allocating huge memory
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param Ulist a list of P covariance matrices for each mixture component
#' @param posterior_weights the JxP posterior probabilities of each mixture component in Ulist for the data
#' @return PosteriorMean JxR matrix of posterior means
#' @return PosteriorSD JxR matrix of posterior (marginal) standard deviations
#' @return NegativeProb JxR matrix of posterior (marginal) probability of being negative
#' @return ZeroProb JxR matrix of posterior (marginal) probability of being zero
#' @return lfsr JxR matrix of local false sign rates
#' @importFrom ashr compute_lfsr
#' @export
compute_posterior_matrices=function(data,Ulist,posterior_weights){
  R=n_conditions(data)
  J=n_effects(data)
  P=length(Ulist)

  # allocate arrays for returned results
  res_post_mean=array(NA,dim=c(J,R))
  res_post_mean2 = array(NA,dim=c(J,R)) #mean squared value
  res_post_zero=array(NA,dim=c(J,R))
  res_post_neg=array(NA,dim=c(J,R))

  # allocate arrays for temporary calculations
  post_mean=array(NA,dim=c(P,R))
  post_mean2 = array(NA,dim=c(P,R)) #mean squared value
  post_zero=array(NA,dim=c(P,R))
  post_neg=array(NA,dim=c(P,R))

  common_cov = is_common_cov(data) # check whether data have common covariance

  if(common_cov){
    V = get_cov(data,1)
    Vinv <- solve(V)
    U1 = lapply(Ulist,
                function(U){posterior_cov(Vinv, U)}) # compute all the posterior covariances
  }

  for(j in 1:J){
    bhat=as.vector(t(data$Bhat[j,]))##turn i into a R x 1 vector
    if(!common_cov){
      V=get_cov(data,j)
      Vinv <- solve(V)
      U1 = lapply(Ulist, function(U){posterior_cov(Vinv, U)}) # compute all the posterior covariances
    }
    for(p in 1:P){
      mu1 <- as.array(posterior_mean(bhat, Vinv, U1[[p]]))
      post_mean[p,]= mu1
      post_mean2[p,] = mu1^2 + diag(U1[[p]]) #diag(U1) is the posterior variance
      post_var = diag(U1[[p]])
      post_neg[p,] = ifelse(post_var==0,0,pnorm(0,mean=mu1,sqrt(diag(U1[[p]])),lower.tail=T))
      post_zero[p,] = ifelse(post_var==0,1,0)
    }
    res_post_mean[j,] = posterior_weights[j,] %*% post_mean
    res_post_mean2[j,] = posterior_weights[j,] %*% post_mean2
    res_post_zero[j,] = posterior_weights[j,] %*% post_zero
    res_post_neg[j,] = posterior_weights[j,] %*% post_neg
  }
  res_post_sd = sqrt(res_post_mean2 - res_post_mean^2)
  res_lfsr = compute_lfsr(res_post_neg,res_post_zero)
  return(list(PosteriorMean = res_post_mean,
              PosteriorSD = res_post_sd,
              lfdr = res_post_zero,
              NegativeProb = res_post_neg,
              lfsr = res_lfsr))
}
