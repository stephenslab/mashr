#' @title Compute posterior matrices (general version)
#' @description Computes posterior matrices without allocating huge memory
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param A the linear transformation matrix, KxR matrix
#' @param Ulist a list of P covariance matrices for each mixture component
#' @param posterior_weights the JxP posterior probabilities of each mixture component in Ulist for the data
#' @return PosteriorMean JxK matrix of posterior means
#' @return PosteriorSD JxK matrix of posterior (marginal) standard deviations
#' @return NegativeProb JxK matrix of posterior (marginal) probability of being negative
#' @return ZeroProb JxK matrix of posterior (marginal) probability of being zero
#' @return lfsr JxK matrix of local false sign rates
#' @importFrom ashr compute_lfsr
#' @importFrom stats pnorm
compute_posterior_matrices_general_R=function(data,A,Ulist,posterior_weights){
  R=n_conditions(data)
  J=n_effects(data)
  P=length(Ulist)
  RA = nrow(A)

  # allocate arrays for returned results
  res_post_mean=array(NA,dim=c(J,RA))
  res_post_mean2 = array(NA,dim=c(J,RA)) #mean squared value
  res_post_zero=array(NA,dim=c(J,RA))
  res_post_neg=array(NA,dim=c(J,RA))

  # allocate arrays for temporary calculations
  post_mean=array(NA,dim=c(P,RA))
  post_mean2 = array(NA,dim=c(P,RA)) #mean squared value
  post_zero=array(NA,dim=c(P,RA))
  post_neg=array(NA,dim=c(P,RA))

  # check if rows of Shat are same, if so,
  # the covariances are same
  common_cov = is_common_cov_Shat(data)

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
      # Transformation for mu
      muA <- A %*% (mu1 * data$Shat_alpha[j,])

      # Transformation for Cov
      covU = data$Shat_alpha[j,] * t(data$Shat_alpha[j,] * U1[[p]])
      pvar = A %*% (covU %*% t(A))
      if(any(pvar < 0)){
        pvar[pvar < 0] = 0
      }

      post_var = diag(pvar) # nrow(A) vector posterior variance
      post_mean[p,]= muA
      post_mean2[p,] = muA^2 + post_var #post_var is the posterior variance

      post_neg[p,] = ifelse(post_var==0,0,pnorm(0,mean=muA,sqrt(post_var),lower.tail=T))
      post_zero[p,] = ifelse(post_var==0,1,0)
    }
    res_post_mean[j,] = posterior_weights[j,] %*% post_mean
    res_post_mean2[j,] = posterior_weights[j,] %*% post_mean2
    res_post_zero[j,] = posterior_weights[j,] %*% post_zero
    res_post_neg[j,] = posterior_weights[j,] %*% post_neg
  }
  res_post_sd = sqrt(res_post_mean2 - res_post_mean^2)
  res_lfsr = compute_lfsr(res_post_neg,res_post_zero)

  if(!is.null(row.names(data$Bhat))){
    row.names(res_post_mean) = row.names(data$Bhat)
    row.names(res_post_sd) = row.names(data$Bhat)
    row.names(res_post_zero) = row.names(data$Bhat)
    row.names(res_post_neg) = row.names(data$Bhat)
    row.names(res_lfsr) = row.names(data$Bhat)
  }
  if(!is.null(row.names(A))){
    colnames(res_post_mean) = row.names(A)
    colnames(res_post_sd) = row.names(A)
    colnames(res_post_zero) = row.names(A)
    colnames(res_post_neg) = row.names(A)
    colnames(res_lfsr) = row.names(A)
  }
  return(list(PosteriorMean = res_post_mean,
              PosteriorSD = res_post_sd,
              lfdr = res_post_zero,
              NegativeProb = res_post_neg,
              lfsr = res_lfsr))
}
