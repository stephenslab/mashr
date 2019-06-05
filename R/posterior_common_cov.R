#' @title Compute posterior matrices (when error covariance V_j is equal for all observations j)
#' @description Computes posterior matrices without allocating huge memory
#' @param data a mash data object, eg as created by \code{mash_set_data} or \code{mash_set_data_contrast}
#' @param A the linear transformation matrix, Q x R matrix. This is used to compute the posterior for Ab.
#' @param Ulist a list of P covariance matrices for each mixture component
#' @param posterior_weights the JxP posterior probabilities of each mixture component in Ulist for the data
#' @param output_posterior_cov whether or not to output posterior covariance matrices for all effects
#' @param posterior_samples the number of samples to be drawn from the posterior distribution of each effect.
#' @param seed a random number seed to use when sampling from the posteriors. It is used when \code{posterior_samples > 0}.
#' @return PosteriorMean JxQ matrix of posterior means
#' @return PosteriorSD JxQ matrix of posterior (marginal) standard deviations
#' @return NegativeProb JxQ matrix of posterior (marginal) probability of being negative
#' @return ZeroProb JxQ matrix of posterior (marginal) probability of being zero
#' @return lfsr JxQ matrix of local false sign rates
#' @return PosteriorCov QxQxJ array of posterior covariance matrices, if the \code{output_posterior_cov = TRUE}
#' @return PosteriorSamples JxQxM array of samples, if the \code{posterior_samples = M > 0}
#' @importFrom ashr compute_lfsr
#' @importFrom stats pnorm rmultinom
#' @importFrom plyr aaply
#' @importFrom mvtnorm rmvnorm
#' @importFrom abind abind
compute_posterior_matrices_common_cov_R=function(data,A, Ulist, posterior_weights, output_posterior_cov = FALSE,
                                                 posterior_samples = 0, seed = 123){
  R = n_conditions(data)
  J = n_effects(data)
  P = length(Ulist)
  Q = nrow(A)

  # allocate arrays for returned results
  res_post_mean=array(0,dim=c(J,Q))
  res_post_mean2 = array(0,dim=c(J,Q)) #mean squared value
  res_post_zero=array(0,dim=c(J,Q))
  res_post_neg=array(0,dim=c(J,Q))

  if(output_posterior_cov){
    # record sum_{p} pi_{jp}(U1[p] + mu_{jp}mu_{jp}^{T})
    post_sec_w_sum = array(0,dim=c(Q, Q, J))
  }

  if(posterior_samples > 0){
    set.seed(seed)
    res_post_samples = vector('list', J)
    Z = apply(posterior_weights, 1, function(p) rowSums(rmultinom(posterior_samples, 1, p))) # P x J matrix
  }

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
      muA <- (mu1 * data$Shat_alpha) %*% t(A) # J by Q matrix
      # Transformation for Cov
      covU = data$Shat_alpha[1,] * t(data$Shat_alpha[1,] * U1[[p]])
      pvar = A %*% (covU %*% t(A))

      post_var = pmax(0,diag(pvar)) # length Q vector posterior variance

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
      if(output_posterior_cov){
        muA_s <- array(0, dim=c(Q,Q,J))
        muA_s[] <- apply(muA, 1, tcrossprod)
        post_sec <- array(0, dim=c(Q,Q,J))
        post_sec[] <- apply(muA_s, 3, '+', pvar)
        post_sec_w_sum <-post_sec_w_sum + aaply(post_sec, c(1,2), function(x,y){x*y}, posterior_weights[,p])
      }

      if(posterior_samples > 0){
        samples_p = lapply(1:J, function(j) if(Z[p,j] > 0){
          rmvnorm(n=Z[p, j], mean = muA[j,], sigma = pvar)
        }else{matrix(0,0,Q)})

        res_post_samples = Map(rbind, res_post_samples, samples_p)
      }
  }

  res_post_sd = sqrt(res_post_mean2 - res_post_mean^2)
  res_lfsr = compute_lfsr(res_post_neg,res_post_zero)

  posterior_matrices = list(PosteriorMean = res_post_mean,
                            PosteriorSD   = res_post_sd,
                            lfdr          = res_post_zero,
                            NegativeProb  = res_post_neg,
                            lfsr          = res_lfsr)

  if(output_posterior_cov){
    muAw_s <- array(0, dim=c(Q,Q,J))
    muAw_s[] <- apply(res_post_mean, 1, tcrossprod)
    res_post_cov = post_sec_w_sum - muAw_s
    posterior_matrices$PosteriorCov = res_post_cov
  }

  if(posterior_samples > 0){
    # shuffle samples
    res_post_samples = lapply(res_post_samples, function(l) l[sample(posterior_samples),])
    res_post_samples = abind(res_post_samples, along = 0) # dim J x M x Q; dim is J x M if Q = 1
    res_post_samples = array(res_post_samples, dim=c(J,posterior_samples,Q))
    res_post_samples = aperm(res_post_samples, c(1,3,2)) # dim J x Q x M
    dimnames(res_post_samples) <- list(rownames(data$Bhat), row.names(A), paste0("sample_",(1:posterior_samples)))
    posterior_matrices$PosteriorSamples = res_post_samples
  }
  return(posterior_matrices)
}
