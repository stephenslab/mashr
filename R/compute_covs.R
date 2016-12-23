#' Compute a list of covariance matrices using supplied functions, and put them all together
#' @param Bhat an N by R matrix of data to be passed to functions fns
#' @param fns a vector of functions to be applied to Bhat to compute covariance matrices
#' @param Ulist optionally, an existing list of matrices to which the newly computed matrices will be added
#' @return a list of covariance matrices
#' export
compute_Ulist = function(Bhat,
                        fns=c(compute_covs_identity,compute_covs_singletons,compute_covs_allones),
                        Ulist = NULL){
  c(Ulist,unlist(lapply(fns,function(f){f(Bhat)}),recursive=FALSE))
}


#' Compute an identity matrix
#' @param Bhat an N by R matrix of data
#' @return a list with one entry, an R by R identity matrix
#' export
compute_covs_identity = function(Bhat){
  R = ncol(Bhat)
  return(list(diag(R)))
}

#' Compute all R singleton matrices corresponding to condition-specific effects
#' @param Bhat an N by R matrix of data
#' @return a list with R entries, the rth entry being an R by R matrix with all 0s except the (r,r) element is 1
#' export
compute_covs_singletons = function(Bhat){
  R = ncol(Bhat)
  nullmatrix = matrix(0,nrow=R,ncol=R)
  Ulist = list()

  for(r in 1:R){
    Ulist[[r]] = nullmatrix
    Ulist[[r]][r,r] = 1
  }
  return(Ulist)
}

#' Compute an R by R matrix of all 1s
#' @param Bhat an N by R matrix of data
#' @return a list with 1 entry, the R by R matrix of all 1s
#' export
compute_covs_allones = function(Bhat){
  R = ncol(Bhat)
  onematrix = matrix(1,nrow=R,ncol=R)
  return(list(onematrix))
}
