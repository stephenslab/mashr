#' @title Compute posterior matrices (when error covariance V_j is equal for all observations j)
#' @description Computes posterior matrices without allocating huge memory
#' @param data a mash data object, eg as created by \code{mash_set_data}
#' @param A the linear transformation matrix, K x R matrix. This is used to compute the posterior for Ab.
#' @param Ulist a list of P covariance matrices for each mixture component
#' @param posterior_weights the JxP posterior probabilities of each mixture component in Ulist for the data
#' @return PosteriorMean JxK matrix of posterior means
#' @return PosteriorSD JxK matrix of posterior (marginal) standard deviations
#' @return NegativeProb JxK matrix of posterior (marginal) probability of being negative
#' @return ZeroProb JxK matrix of posterior (marginal) probability of being zero
#' @return lfsr JxK matrix of local false sign rates
#' @importFrom ashr compute_lfsr
#' @importFrom stats pnorm
compute_posterior_matrices_common_cov_R=function(data,A, Ulist, posterior_weights){
  R = n_conditions(data)
  J = n_effects(data)
  P = length(Ulist)
  K = nrow(A)

  # allocate arrays for returned results
  res_post_mean=array(0,dim=c(J,K))
  res_post_mean2 = array(0,dim=c(J,K)) #mean squared value
  res_post_zero=array(0,dim=c(J,K))
  res_post_neg=array(0,dim=c(J,K))

  if((!is_common_cov_Shat(data)) && (!is_common_cov_Shat_alpha(data))){
    stop("Can't call common_cov routine with data where covariance varies")
  } # check whether data have common covariance

  V = get_cov(data,1)
  Vinv <- solve(V)
  # compute all the posterior covariances
  U1 = lapply(Ulist,function(U){posterior_cov(Vinv, U)})

  for(p in 1:P){
      mu1 <- posterior_mean_matrix(data$Bhat, Vinv, U1[[p]]) # J by R matrix
      # Transformation for mu
      muA <- (mu1 * data$Shat_alpha) %*% t(A) # J by nrow(A) matrix
      # Transformation for Cov
      covU = data$Shat_alpha[1,] * t(data$Shat_alpha[1,] * U1[[p]])
      pvar = A %*% (covU %*% t(A))
      # correct for computing error
      if(any(pvar < 0)){
        pvar[pvar < 0] = 0
      }
      post_var = diag(pvar) # nrow(A) vector posterior variance

      res_post_mean= res_post_mean + posterior_weights[,p] * muA

      res_post_mean2 = res_post_mean2 +
          posterior_weights[,p] * t(t(muA^2) + post_var)

      null.cond = (post_var==0) # which conditions are null
      res_post_zero[,null.cond] = res_post_zero[,null.cond] +
        posterior_weights[,p]

      # only the non-null conditions have a probability of being negative
      if (!all(null.cond))
        res_post_neg[,!null.cond] = res_post_neg[,!null.cond] +
          posterior_weights[,p] * t(
            pnorm(0,mean=t(muA[,!null.cond]),sqrt(post_var[!null.cond]),
                  lower.tail=T))

  }

  res_post_sd = sqrt(res_post_mean2 - res_post_mean^2)
  res_lfsr = compute_lfsr(res_post_neg,res_post_zero)

  # Add names
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
