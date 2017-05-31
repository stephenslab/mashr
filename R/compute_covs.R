#' Compute a list of canonical covariance matrices
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param cov_methods a vector of strings indicating the matrices to be used:
#' "identity" for the identity (effects are independent among conditions);
#' "singletons" for the set of matrices with just one non-zero entry x_{jj} = 1 (j=1,...,R); (effect specific to condition j);
#' "equal_effects" for the matrix of all 1s (effects are equal among conditions);
#' "simple_het" for a set of matrices with 1s on the diagonal and all off-diagonal elements equal to 0.25, 0.5 or 0.75; see \code{cov_simple_het} for details; (effects are correlated among conditions).
#' @return a list of covariance matrices
#' @details The default is that this function computes covariance matrices corresponding
#' to the "bmalite" models.
#' @examples data = set_mash_data(Bhat = cbind(c(1,2),c(3,4)), Shat = cbind(c(1,1),c(1,1)))
#'  cov_canonical(data)
#'  cov_canonical(data,"singletons")
#'  cov_canonical(data,c("id","sing")) # can use partial matching of names
#' @importFrom utils modifyList
#' @export
cov_canonical = function(data,
                       cov_methods= c("identity","singletons","equal_effects","simple_het")){
  res = list()
  for(i in 1:length(cov_methods)){
    if(is.character(cov_methods[[i]])){
      name = match.arg(cov_methods[[i]], names(cov_methods()))
      cov_method = cov_methods()[[ name ]]
    } else {
      cov_method = cov_methods[[i]]
      name = names(cov_methods)[i]
    }
    res[[i]] =  do.call(cov_method$fn, args = modifyList(list(data=data), as.list(cov_method$args) ))

    if(is.list(res[[i]])){ #append _i to names if function returns multiple matrices
      names(res[[i]]) = paste0(name, "_", 1:length(res[[i]]) )
    } else {
      res[[i]] = list(res[[i]])
      names(res[[i]]) = name
    }
  }
  res = unlist(res,recursive=FALSE)
  return(res)
}

# #' Compute a list of data-driven covariance matrices;
# #' @param data a mash data object, eg as created by \code{set_mash_data}
# #' @return a list of covariance matrices
# #' @details currently just applies Extreme Deconvolution to a subset of signals; future expansion expected
# #' @export
# cov_data_driven = function(data,subset=NULL){
#   data2cov_ed(data,subset=subset)
# }


#' A function that maps names of covariance matrices to covariance computation functions
#' @return a named list of defaults for covariance function calculation
cov_methods = function(){
  list("null"= list(fn = cov_all_zeros, args = NULL),
        "identity" = list(fn = cov_identity, args= NULL),
        "singletons" = list(fn = cov_singletons, args=NULL),
        "equal_effects" = list(fn = cov_equal_effects, args = NULL),
        "simple_het" = list(fn = cov_simple_het, args = list(corr=c(0.25,0.5,0.75)))
  )
}


#' Compute an identity matrix
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return the R by R identity matrix
cov_identity = function(data){
  R = n_conditions(data)
  diag(R)
}

#' Compute all R singleton matrices corresponding to condition-specific effects
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with R entries, the rth entry being an R by R matrix with all 0s except the (r,r) element is 1
cov_singletons = function(data){
  R = n_conditions(data)
  nullmatrix = matrix(0,nrow=R,ncol=R)
  Ulist = list()

  for(r in 1:R){
    Ulist[[r]] = nullmatrix
    Ulist[[r]][r,r] = 1
  }
  return(Ulist)
}

#' Compute all the singleton matrices corresponding to condition-specific effects in first condition only; used for testing purposes
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return an R by R matrix with all 0s except the (1,1) element is 1
cov_first_singleton = function(data){
  R = n_conditions(data)
  res = matrix(0,nrow=R,ncol=R)
  res[1,1]=1
  return(res)
}

#' Compute an R by R matrix of all 1s
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with 1 entry, the R by R matrix of all 1s
cov_equal_effects = function(data){
  R = n_conditions(data)
  matrix(1,nrow=R,ncol=R)
}

#' Compute an R by R matrix of all 0s
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with 1 entry, the R by R matrix of all 0s
cov_all_zeros = function(data){
  R = ncol(data$Bhat)
  matrix(0,nrow=R,ncol=R)
}

#' Compute covariance matrices with diagonal element 1 and off-diagonal element corr
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param corr a vector containing the correlations to be used
#' @return a list of matrices, one for each value in corr
#' @export
cov_simple_het = function(data, corr=c(0.25,0.5,0.75)){
  R = n_conditions(data)
  simplehet=list()
  for(i in 1:length(corr)){
    if(corr[i]>1 | corr[i]<(-1)){stop("corr must be between -1 and 1")}
    simplehet[[i]] = matrix(corr[i],nrow=R,ncol=R)
    diag(simplehet[[i]]) <- 1
  }
  return(simplehet)
}

#' Scale each covariance matrix in list Ulist by a scalar in vector grid
#' @param Ulist a list of matrices
#' @param grid a vector of scaling factors (standard deviaions)
#' @return a list with length length(Ulist)*length(grid), with values grid[i]^2*Ulist[[j]]
#' @export
scale_cov = function(Ulist, grid){
  orig_names = names(Ulist)
  Ulist = unlist( lapply(grid^2, function(x){multiply_list(Ulist,x)}), recursive=FALSE)
  names(Ulist) = unlist( lapply(1:length(grid), function(x){paste0(orig_names,".",x)}), recursive=FALSE)
  return(Ulist)
}

#' multiply each element of a list by scalar
#' (In our application each element of the list is a matrix)
multiply_list = function(Ulist, x){lapply(Ulist, function(U){x*U})}


#' normalize a covariance matrix soits maximum diagonal element is 1
#' Divides each element of the matrix by its maximum diagonal element (provided that it is non-zero)
normalize_cov = function(U){
  if(max(diag(U))!=0){
    U = U/max(diag(U))
  }
  return(U)
}

#' normalize a list of covariance matrices by
#' applying \link{\code{normalize_cov}} to each element
normalize_Ulist = function(Ulist){lapply(Ulist,normalize_cov)}

#' Create names for covariance matrices
#' @param names a string
#' @param suffixes
#' Adds _suffixes to names for each element of suffixes
make_names = function(names,suffixes){paste0(names,"_",suffixes)}

check_dim = function(mat,R){
  if(!identical(dim(mat),c(R,R))){stop("Dimension of matrix must be R by R")}
}

