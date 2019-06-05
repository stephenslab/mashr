#' @title Compute posterior matrices (general version)
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
compute_posterior_matrices_general_R=function(data,A,Ulist,posterior_weights,output_posterior_cov = FALSE,
                                              posterior_samples = 0, seed = 123){
  R=n_conditions(data)
  J=n_effects(data)
  P=length(Ulist)
  Q = nrow(A)

  # allocate arrays for returned results
  res_post_mean=array(NA,dim=c(J,Q))
  res_post_mean2 = array(NA,dim=c(J,Q)) #mean squared value
  res_post_zero=array(NA,dim=c(J,Q))
  res_post_neg=array(NA,dim=c(J,Q))

  # allocate arrays for temporary calculations
  post_mean=array(NA,dim=c(P,Q))
  post_mean2 = array(NA,dim=c(P,Q)) #mean squared value
  post_zero=array(NA,dim=c(P,Q))
  post_neg=array(NA,dim=c(P,Q))

  if(output_posterior_cov){
    res_post_cov = array(NA, dim=c(Q, Q, J))
    post_cov=array(NA, dim=c(Q,Q,P))
  }

  if(posterior_samples > 0){
    set.seed(seed)
    # length J list
    res_post_samples = vector("list", J)
  }

  # check if rows of Shat are same, if so,
  # the covariances are same
  common_cov = is_common_cov_Shat(data)

  if(data$commonV && common_cov){
    V = get_cov(data,1)
    Vinv <- solve(V)
    U1 = lapply(Ulist,
                function(U){posterior_cov(Vinv, U)}) # compute all the posterior covariances
  }

  for(j in 1:J){
    bhat=as.vector(t(data$Bhat[j,]))##turn j into a R x 1 vector
    if(!common_cov || !data$commonV){
      V=get_cov(data,j)
      Vinv <- solve(V)
      U1 = lapply(Ulist, function(U){posterior_cov(Vinv, U)}) # compute all the posterior covariances
    }
    if(posterior_samples > 0){
      samples_j = vector('list', P)
      z = rowSums(rmultinom(posterior_samples, 1, posterior_weights[j,]))
    }
    for(p in 1:P){
      mu1 <- as.array(posterior_mean(bhat, Vinv, U1[[p]]))
      # Transformation for mu
      muA <- A %*% (mu1 * data$Shat_alpha[j,])

      # Transformation for Cov
      covU = data$Shat_alpha[j,] * t(data$Shat_alpha[j,] * U1[[p]])
      pvar = A %*% (covU %*% t(A))

      post_var = pmax(0,diag(pvar)) # Q vector posterior variance
      post_mean[p,]= muA
      post_mean2[p,] = muA^2 + post_var #post_var is the posterior variance

      post_neg[p,] = ifelse(post_var==0,0,pnorm(0,mean=muA,sqrt(post_var),lower.tail=T))
      post_zero[p,] = ifelse(post_var==0,1,0)

      if(output_posterior_cov){
        post_cov[,,p] = pvar + tcrossprod(muA)
      }

      if(posterior_samples > 0){
        if(z[p] > 0){
          samples_j[[p]] = rmvnorm(z[p], mean=muA, sigma = pvar)
        } else{
          samples_j[[p]] = matrix(0,0,Q)
        }
      }
    }
    res_post_mean[j,] = posterior_weights[j,] %*% post_mean
    res_post_mean2[j,] = posterior_weights[j,] %*% post_mean2
    res_post_zero[j,] = posterior_weights[j,] %*% post_zero
    res_post_neg[j,] = posterior_weights[j,] %*% post_neg
    if(output_posterior_cov){
      res_post_cov[,,j] = aaply(post_cov, c(1,2), function(x,y){crossprod(x,y)}, posterior_weights[j,]) - tcrossprod(res_post_mean[j,])
    }
    if(posterior_samples > 0){
      samples_j_full = do.call(rbind, samples_j)
      res_post_samples[[j]] = samples_j_full[sample(posterior_samples),]
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
    posterior_matrices$PosteriorCov = res_post_cov
  }
  if(posterior_samples > 0){
    res_post_samples = abind(res_post_samples, along = 0,force.array=TRUE) # dim J x M x Q; dim is J x M if Q = 1
    res_post_samples = array(res_post_samples, dim=c(J,posterior_samples,Q))
    res_post_samples = aperm(res_post_samples, c(1,3,2)) # dim J x Q x M
    dimnames(res_post_samples) <- list(rownames(data$Bhat), row.names(A), paste0("sample_",(1:posterior_samples)))

    posterior_matrices$PosteriorSamples = res_post_samples
  }
  return(posterior_matrices)
}
