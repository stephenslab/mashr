#' @title Compute posterior matrices (when error covariance V_j is equal for all observations j)
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
#' @importFrom stats pnorm
#' @export
compute_posterior_matrices_common_cov_R=function(data,Ulist,posterior_weights){
  R=n_conditions(data)
  J=n_effects(data)
  P=length(Ulist)

  # allocate arrays for returned results
  res_post_mean=array(0,dim=c(J,R))
  res_post_mean2 = array(0,dim=c(J,R)) #mean squared value
  res_post_zero=array(0,dim=c(J,R))
  res_post_neg=array(0,dim=c(J,R))

  if(!is_common_cov(data)){
    stop("Can't call common_cov routine with data where covariance varies")
  } # check whether data have common covariance

  V = get_cov(data,1)
  Vinv <- solve(V)
  # compute all the posterior covariances
  U1 = lapply(Ulist,function(U){posterior_cov(Vinv, U)})

  for(p in 1:P){
      mu1 <- posterior_mean_matrix(data$Bhat, Vinv, U1[[p]]) # J by R matrix
      post_var = diag(U1[[p]]) # R vector

      res_post_mean= res_post_mean + posterior_weights[,p] * mu1

      res_post_mean2 = res_post_mean2 +
          posterior_weights[,p] * t(t(mu1^2) + post_var)

      null.cond = (post_var==0) # which conditions are null
      res_post_zero[,null.cond] = res_post_zero[,null.cond] +
        posterior_weights[,p]

      # only the non-null conditions have a probability of being negative
      res_post_neg[,!null.cond] = res_post_neg[,!null.cond] +
        posterior_weights[,p] * t(
          pnorm(0,mean=t(mu1[,!null.cond]),sqrt(post_var[!null.cond]),
                lower.tail=T))

  }

  res_post_sd = sqrt(res_post_mean2 - res_post_mean^2)
  res_lfsr = compute_lfsr(res_post_neg,res_post_zero)
  return(list(PosteriorMean = res_post_mean,
              PosteriorSD = res_post_sd,
              lfdr = res_post_zero,
              NegativeProb = res_post_neg,
              lfsr = res_lfsr))
}
