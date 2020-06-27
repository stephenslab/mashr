#' @title Perform PCA on data and return list of candidate covariance
#' matrices
#'
#' @param data a mash data object
#'
#' @param npc the number of PCs to use
#'
#' @param subset indices of the subset of data to use (set to NULL for
#' all data)
#'
#' @return Returns a list of covariance matrices: the npc rank-one
#' covariance matrices based on the first npc PCs, and the rank npc
#' covariance matrix
#' @examples
#' data = mash_set_data(Bhat = cbind(c(1,2),c(3,4)), Shat = cbind(c(1,1),c(1,1)))
#' cov_pca(data,2)
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
cov_pca = function(data,npc,subset = NULL){
  assert_that(npc>1)
  assert_that(npc<=n_conditions(data))
  if(is.null(subset)){subset = 1:n_effects(data)}
  res.svd = svd(data$Bhat[subset,],nv=npc,nu=npc)
  # FIXME: we need to think of for the EE case what to use for svd input: Bhat or Bhat/Shat
  f = res.svd$v
  Ulist = cov_from_factors(t(f), "PCA")
  d = diag(res.svd$d[1:npc])
  Ulist = c(Ulist, list("tPCA"= f %*% d^2 %*% t(f)/length(subset)))
  return(Ulist)
}

#' @title Perform "extreme deconvolution" (Bovy et al) on a subset of
#' the data
#'
#' @param data a mash data object
#'
#' @param Ulist_init a named list of covariance matrices to use to
#' initialize ED; default is to use matrices from PCs
#'
#' @param subset a subset of data to be used when ED is run (set to
#' NULL for all the data)
#'
#' @param algorithm algorithm to run ED
#'
#' @param ... other arguments to be passed to ED algorith, see
#' \code{\link{extreme_deconvolution}} for algorithm 'bovy', or
#' \code{\link{teem_wrapper}} for algorithm 'teem'
#'
#' @details Runs the extreme deconvolution algorithm from Bovy et al
#' (Annals of Applied Statistics) to estimate data-driven covariance
#' matrices. It can be initialized with, for example running \code{cov_pca} with, 
#' say, 5 PCs.
#' @examples
#' data = mash_set_data(Bhat = cbind(c(1,2),c(3,4)), Shat = cbind(c(1,1),c(1,1)))
#' U_pca = cov_pca(data,2)
#' U_x = apply(data$Bhat, 2, function(x) x - mean(x))
#' U_xx = t(U_x) %*% U_x / nrow(U_x)
#' cov_ed(data,c(U_pca, list(xx = U_xx)))
#'
#' @export
#'
cov_ed = function(data, Ulist_init, subset = NULL,
                  algorithm=c('bovy', 'teem'), ...) {
  algorithm = match.arg(algorithm)
  if (algorithm=='bovy') {
    Ulist_ed = bovy_wrapper(data, Ulist_init, subset, ...)$Ulist
  } else {
    Ulist_ed = teem_wrapper(data, Ulist_init, subset, ...)$U
  }
  names(Ulist_ed) = make_names("ED", if(is.null(names(Ulist_init))) 1:length(Ulist_ed) else names(Ulist_init))
  Ulist_ed
}

# For a vector x, return the rank one matrix xx'.
r1cov=function(x){x %*% t(x)}

#' produce list of rank-1 covariance matrices corresponding to rows of f
#'
#' @param f a matrix of factors (each row is a factor)
#'
#' @param name a string indicating the name to use
#'
#' @return a list of rank one matrices whose kth element is f[k,]
#' f[k,]' and named name_k
#'
#' @keywords internal
#'
cov_from_factors = function(f, name){
  Ulist = list()
  for(i in 1:nrow(f)){
    Ulist = c(Ulist,list(r1cov(f[i,])))
  }
  names(Ulist) = paste0(name,"_",(1:nrow(f)))
  return(Ulist)
}
