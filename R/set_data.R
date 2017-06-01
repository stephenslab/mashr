#' Create a data object for mash analysis
#' @param Bhat an N by R matrix of observed estimates
#' @param Shat an N by R matrix of corresponding standard errors
#' @param V an R by R correlation matrix of error correlations
#'  [So Bhat_j distributed as N(B_j,diag(Shat_j) V diag(Shat_j)) where _j denotes the jth row of a matrix]
#'  Defaults to identity
#' @return a data object for passing into mash functions
#' @export
set_mash_data = function(Bhat,Shat,V=diag(ncol(Bhat))){
  return(list(Bhat=Bhat, Shat=Shat, V=V))
}

#' return the covariance matrix for jth data point from mash data object
get_cov = function(data,j){
  data$Shat[j,] * t(data$V * data$Shat[j,]) # quicker than diag(Shat[j,]) %*% V %*% diag(Shat[j,])
}

n_conditions = function(data){ncol(data$Bhat)}

n_effects = function(data){nrow(data$Bhat)}
