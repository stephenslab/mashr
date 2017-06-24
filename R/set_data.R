#' Create a data object for mash analysis
#' @param Bhat an N by R matrix of observed estimates
#' @param Shat an N by R matrix of corresponding standard errors. Shat can be a scalar if all standard errors are equal.
#' This is most useful if Bhat is a matrix of Z scores, so elements of Shat are all 1.
#' @param V an R by R correlation matrix of error correlations; must be positive definite.
#'  [So Bhat_j distributed as N(B_j,diag(Shat_j) V diag(Shat_j)) where _j denotes the jth row of a matrix]
#'  Defaults to identity
#' @return a data object for passing into mash functions
#' @export
set_mash_data = function(Bhat,Shat,V=diag(ncol(Bhat))){
  if(length(Shat)==1){Shat = matrix(Shat,nrow=nrow(Bhat),ncol=ncol(Bhat))}
  if(!identical(dim(Bhat),dim(Shat))){stop("dimensions of Bhat and Shat must match")}
  if(det(V)<1e-15){stop("V must be positive definite")}
  return(list(Bhat=Bhat, Shat=Shat, V=V))
}

# Return the covariance matrix for jth data point from mash data object.
get_cov = function(data,j){
  data$Shat[j,] * t(data$V * data$Shat[j,]) # quicker than diag(Shat[j,]) %*% V %*% diag(Shat[j,])
}

n_conditions = function(data){ncol(data$Bhat)}

n_effects = function(data){nrow(data$Bhat)}
