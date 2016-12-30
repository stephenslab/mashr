#' Compute a list of covariance matrices using supplied functions, and put them all together
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param fns a vector of functions to be applied to Bhat to compute covariance matrices
#' @param Ulist optionally, an existing list of matrices to which the newly computed matrices will be added
#' @return a list of covariance matrices
#' @export
compute_Ulist = function(data,
                        fns=c(cov_identity,cov_singletons,cov_allones,cov_allzeros),
                        Ulist = NULL){
  if(length(fns)==1){fns = c(fns)}
  names(fns)=NULL #avoid any names that happen to be passed in being passed on to the output
  c(Ulist,unlist(lapply(fns,function(f){f(data)}),recursive=FALSE))
}

#maps string names to functions used to compute different types of covariance matrix
Umap = function(){
  list("allzeros"= cov_allzeros,
       "id" = cov_identity,
       "singletons" = cov_singletons,
       "allones" = cov_allones)
}

#' Compute a list of covariance matrices (with methods indicated by names)
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param Utypes a vector of strings that indicate which covariance matrices to compute. See \code{Umap}.
#' @param Ulist optionally, an existing list of matrices to which the newly computed matrices will be added
#' @return a list of covariance matrices
#' @export
compute_Ulist_byname = function(data,
           Utypes=c("id","singletons","allones","allzeros"),
           Ulist = NULL){
  Utypes = match.arg(Utypes,several.ok=TRUE)
  Ufns = Umap()[Utypes] # map the names in Utypes to functions
  compute_Ulist(data, Ufns, Ulist)
}

#' Compute an identity matrix
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with one entry, an R by R identity matrix
#' @export
cov_identity = function(data){
  R = n_conditions(data)
  return(list(id=diag(R)))
}

#' Compute all R singleton matrices corresponding to condition-specific effects
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with R entries, the rth entry being an R by R matrix with all 0s except the (r,r) element is 1
#' @export
cov_singletons = function(data){
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
cov_allones = function(data){
  R = n_conditions(data)
  onematrix = matrix(1,nrow=R,ncol=R)
  return(list(allones = onematrix))
}

#' Compute an R by R matrix of all 0s
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with 1 entry, the R by R matrix of all 0s
#' @export
cov_allzeros = function(data){
  R = ncol(data$Bhat)
  zeromatrix = matrix(0,nrow=R,ncol=R)
  return(list(allzeros = zeromatrix))
}

#' Scale each covariance matrix in list Ulist by a scalar in vector grid
#' @param Ulist a list of matrices
#' @param grid a vector of scaling factors
#' @return a list with length length(Ulist) \times length(grid), with values grid[i]\times Ulist[[j]]
#' @export
scale_cov = function(Ulist, grid){
  orig_names = names(Ulist)
  Ulist = unlist( lapply(grid, function(x){multiply_list(Ulist,x)}), recursive=FALSE)
  names(Ulist) = unlist( lapply(1:length(grid), function(x){paste0(orig_names,".",x)}), recursive=FALSE)
  return(Ulist)
}

#' multiply each element of a list by scalar
#' (In our application each element of the list is a matrix)
multiply_list = function(Ulist, x){lapply(Ulist, function(U){x*U})}

