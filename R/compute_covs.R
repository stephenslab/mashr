#' Compute a list of covariance matrices using supplied functions, and put them all together
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @param cov_methods a named list defining the methods (functions and arguments) to be used
#' @param Ulist optionally, an existing list of matrices to which the newly computed matrices will be added
#' @return a list of covariance matrices
#' @details Each element of methods may be either a character string which names a pre-defined method (see \code{cov_methods}) or a list with elements fn and args, specifying the function to call, and any additional arguments
#' @examples data = set_mash_data(Bhat = cbind(c(1,2),c(3,4)), Shat = cbind(c(1,1),c(1,1)))
#'  compute_cov(data,"id")
#'  compute_cov(data,"singletons")
#'  compute_cov(data,c("id","sing")) # can use partial matching of names
#'  compute_cov(data,list("sing","id",eg = list(fn = cov_simple_het, args= list(corr=c(0.1,0.2)))))
#' @export
compute_cov = function(data,
                        cov_methods,
                        Ulist = NULL){
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

    if(is.list(res[[i]])){ #append _i to names if function retuns multiple matrices
      names(res[[i]]) = paste0(name, "_", 1:length(res[[i]]) )
    } else {
      res[[i]] = list(res[[i]])
      names(res[[i]]) = name
    }
  }
  res = unlist(res,recursive=FALSE)
  if(!is.null(Ulist)){res = c(Ulist,res)}
  return(res)
}



#' A function that maps names of covariance matrices to covariance computation functions
#' @return a named list of defaults for covariance function calculation
cov_methods = function(){
  list("null"= list(fn = cov_all_zeros, args = NULL),
        "identity" = list(fn = cov_identity, args= NULL),
        "singletons" = list(fn = cov_singletons, args=NULL),
        "all_ones" = list(fn = cov_all_ones, args = NULL),
        "simple_het" = list(fn = cov_simple_het, args = list(corr=c(0.25,0.5,0.75)))
  )
}


#' Compute an identity matrix
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return the R by R identity matrix
#' @export
cov_identity = function(data){
  R = n_conditions(data)
  diag(R)
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
  return(Ulist)
}

#' Compute all the singleton matrices corresponding to condition-specific effects in first condition only; used for testing purposes
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return an R by R matrix with all 0s except the (1,1) element is 1
#' @export
cov_first_singleton = function(data){
  R = n_conditions(data)
  res = matrix(0,nrow=R,ncol=R)
  res[1,1]=1
  return(res)
}

#' Compute an R by R matrix of all 1s
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with 1 entry, the R by R matrix of all 1s
#' @export
cov_all_ones = function(data){
  R = n_conditions(data)
  matrix(1,nrow=R,ncol=R)
}

#' Compute an R by R matrix of all 0s
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list with 1 entry, the R by R matrix of all 0s
#' @export
cov_all_zeros = function(data){
  R = ncol(data$Bhat)
  matrix(0,nrow=R,ncol=R)
}

#' For each element of corr compute a matrix with diagonal element 1 and off-diagonal element corr
#' @param data a mash data object, eg as created by \code{set_mash_data}
#' @return a list of matrices
#' @export
cov_simple_het = function(data, corr){
  R = n_effects(data)
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
