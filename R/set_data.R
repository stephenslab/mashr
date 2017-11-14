#' Create a data object for mash analysis
#' @param Bhat an N by R matrix of observed estimates
#' @param Shat an N by R matrix of corresponding standard errors. Shat can be a scalar if all standard errors are equal.
#' This is most useful if Bhat is a matrix of Z scores, so elements of Shat are all 1. Default is 1.
#' @param alpha Numeric value of alpha parameter in the model. alpha = 0 for Exchangeable Effects (EE), alpha = 1 for Exchangeable Z-scores (EZ). Default is 1. Please refer to equation (3.2) of M. Stephens 2016, Biostatistics for a discussion on alpha.
#' @param df an N by R matrix of corresponding degrees of freedom of the t-statistic Bhat/Shat. Can be a scalar if all degrees of freedom are equal. Default is inf (for large samples).
#' @param pval an N by R matrix of p-values of t-statistic Bhat/Shat. Shat and df should not be specified when pval is provided.
#' @param V an R by R correlation matrix of error correlations; must be positive definite.
#'  [So Bhat_j distributed as N(B_j,diag(Shat_j) V diag(Shat_j)) where _j denotes the jth row of a matrix]
#'  Defaults to identity
#' @return a data object for passing into mash functions
#' @export
set_mash_data = function(Bhat,Shat=NULL,alpha,df=NULL,pval=NULL,V=diag(ncol(Bhat))){
  if (missing(Shat) && missing(pval)) {
    Shat = 1
  }
  if (!missing(pval) && !missing(Shat)) {
    stop("Either Shat or pval can be specified but not both.")
  }
  if (!missing(pval) && !missing(df)) {
    stop("Either df or pval can be specified but not both.")
  }
  if (!missing(pval)) {
    ## Shat and df have to be NULL
    Shat = Bhat / p2z(pval, Bhat)
  }
  if(length(Shat)==1){Shat = matrix(Shat,nrow=nrow(Bhat),ncol=ncol(Bhat))}
  if(!identical(dim(Bhat),dim(Shat))){stop("dimensions of Bhat and Shat must match")}
  if(det(V)<1e-15){stop("V must be positive definite")}
  if(missing(df) || is.infinite(df)) {
    Shat_orig = NULL
  } else {
    if(length(df)==1){df = matrix(df,nrow=nrow(Bhat),ncol=ncol(Bhat))}
    ## Shat = Bhat/Z where Z is the Z score corresponding to a p value from a t test done on (Bhat,Shat_orig,df)
    Shat_orig = Shat
    Shat = Bhat / p2z(2 * pt(-abs(Bhat/Shat), df), Bhat)
  }

  # transform data according to alpha
  if (alpha != 0 && !all(Shat == 1)) {
    ## alpha models dependence of effect size on standard error
    ## alpha > 0 implies larger effects has large standard error
    ## a special case when alpha = 1 is the EZ model
    Shat_alpha = Shat^alpha
    Bhat = Bhat / Shat_alpha
    Shat = Shat^(1-alpha)
  } else {
    Shat_alpha = matrix(1, nrow(Shat), ncol(Shat))
  }
  return(list(Bhat=Bhat, Shat=Shat, V=V, alpha=alpha))
}

# Return the covariance matrix for jth data point from mash data object.
get_cov = function(data,j){
  data$Shat[j,] * t(data$V * data$Shat[j,]) # quicker than diag(Shat[j,]) %*% V %*% diag(Shat[j,])
}

n_conditions = function(data){ncol(data$Bhat)}

n_effects = function(data){nrow(data$Bhat)}

p2z = function(pval, Bhat) {
  z = abs(qnorm(pval / 2))
  z[which(Bhat < 0)] = -1 * z[which(Bhat < 0)]
  return(z)
}
