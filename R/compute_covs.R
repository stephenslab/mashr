#' Compute a list of covariance matrices using supplied functions, and put them all together
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param fns a vector of functions to be applied to Bhat to compute covariance matrices
#' @param Ulist optionally, an existing list of matrices to which the newly computed matrices will be added
#' @return a list of covariance matrices
#' @export
compute_Ulist = function(data,
                        fns=c(compute_covs_identity,compute_covs_singletons,compute_covs_allones),
                        Ulist = NULL){
  if(length(fns)==1){fns = c(fns)}
  c(Ulist,unlist(lapply(fns,function(f){f(data)}),recursive=FALSE))
}


#' Compute an identity matrix
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with one entry, an R by R identity matrix
#' @export
compute_covs_identity = function(data){
  R = n_conditions(data)
  return(list(id=diag(R)))
}

#' Compute all R singleton matrices corresponding to condition-specific effects
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with R entries, the rth entry being an R by R matrix with all 0s except the (r,r) element is 1
#' @export
compute_covs_singletons = function(data){
  R = n_conditions(data)
  nullmatrix = matrix(0,nrow=R,ncol=R)
  Ulist = list()

  for(r in 1:R){
    Ulist[[r]] = nullmatrix
    Ulist[[r]][r,r] = 1
  }
  names(Ulist) <- paste0("singleton_", 1:R)
  return(Ulist)
}

#' Compute an R by R matrix of all 1s
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with 1 entry, the R by R matrix of all 1s
#' @export
compute_covs_allones = function(data){
  R = n_conditions(data)
  onematrix = matrix(1,nrow=R,ncol=R)
  return(list(onematrix = onematrix))
}

#' Compute an R by R matrix of all 0s
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with 1 entry, the R by R matrix of all 0s
#' @export
compute_covs_allzeros = function(data){
  R = ncol(data$Bhat)
  zeromatrix = matrix(0,nrow=R,ncol=R)
  return(list(zeromatrix = zeromatrix))
}

#' Scale each matrix in list Ulist by a scalar in vector grid
#' @param Ulist a list of matrices
#' @param grid a vector of scaling factors
#' @return a list with length length(Ulist) \times length(grid), with values grid[i]\times Ulist[[j]]
#' @export
scale_Ulist = function(Ulist, grid){
  orig_names = names(Ulist)
  Ulist = unlist( lapply(grid, function(x){multiply_list(Ulist,x)}), recursive=FALSE)
  names(Ulist) = unlist( lapply(1:length(grid), function(x){paste0(orig_names,".",x)}), recursive=FALSE)
  return(Ulist)
}

#' multiply each element of a list by scalar
#' (In our application each element of the list is a matrix)
multiply_list = function(Ulist, x){lapply(Ulist, function(U){x*U})}

