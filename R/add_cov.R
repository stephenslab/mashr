#' Add a list of covariance matrices to m
#' @param m a mash object
#' @param Ulist a list of covariance matrices to add to m
#' @export
mash_add_cov = function(m, cov_methods){
  m$Ulist = compute_cov(m$data,cov_methods,m$Ulist)
}


#' Add singleton covariance matrices to m
#' @param m a mash object
#' @export
mash_add_cov_singletons = function(m){mash_add_cov(m,"singletons")}

#' Add identify covariance matrix to m
#' @param m a mash object
#' @export
mash_add_cov_identity = function(m){mash_add_cov(m,"identity")}

#' Add null covariance matrix to m
#' @param m a mash object
#' @export
mash_add_cov_null = function(m){mash_add_cov(m,"null")}

#' Add a covariance matrix of equal effects to m
#' @param m a mash object
#' @export
mash_add_cov_equal_effects = function(m){mash_add_cov(m,"equal_effects")}

#' Add covariance matrices allowing for heterogenous (positively correlated) effects
#' @param m a mash object
#' @details These covariance matrices assume all pairs of conditions are equally correlated.
#' @export
mash_add_cov_simple_het = function(m){mash_add_cov(m,"simple_het")}


#' For a vector x, return the rank one matrix xx'
r1cov=function(x){x %*% t(x)}

#' Add rank 1 covariance matrices to m
#' @param m a mash object
#' @param f a matrix, whose rows are each of length R; if f is a vector, co-erced to a matrix with one row
#' @details For each row x of matrix f, this adds the rank one matrix xx' as a covariance matrix to m
#' @export
mash_add_cov_r1=function(m,f){
  if(is.vector(f)){f = matrix(f,nrow=1)}
  nr1 = length(grep("rank1", list_cov(m))) # find out how many rank 1 matrices currently exist
  to_add = list()
  for(i in 1:nrow(f)){
    to_add = c(to_add,list(r1cov(f[i,])))
  }
  names(to_add) = paste0("rank1_",(1:nrow(f))+nr1)
  mash_add_cov_list(m, to_add)
}

#' Add covariance matrices by doing PCA on data in mash object
#' @param m a mash object
#' @param k number of PCs to use
#' @param subset subset of data points to use in PCA (default of NULL causes top hits in m to be used)
#' @export
mash_add_cov_pca=function(m, k=5, subset=NULL){
  m.svd = svd(get_Z(m,subset))
  f = m.svd$v[,1:k]
  mash_add_cov_r1(m,f)

}


#' Add a list of covariance matrices to m
#' @param m a mash object
#' @param Ulist a named list of covariance matrices to add to m
#' @details This allows the user to provide their own list of covariance matrices if they desire,
#' in addition to built-in covariance options
#' @export
mash_add_cov_list = function(m, Ulist){
  if(!is.list(Ulist) | is.null(names(Ulist))){stop("Ulist must be a named list")}
  lapply(Ulist, check_dim, R=n_conditions.mash(m))
  m$Ulist = c(m$Ulist, as.list(Ulist))
}


check_dim = function(mat,R){
  if(!identical(dim(mat),c(R,R))){stop("Dimension of matrix must be R by R")}
}


#' Return some Bhat values from m
#' @param m a mash object
#' @export
get_Bhat = function(m, subset=NULL){
  if(is.null(subset)){subset = m$tophits}
  return(m$Bhat[subset,])
}

#' Return matrix of Z scores from from m
#' @param m a mash object
#' @export
get_Z = function(m, subset=NULL){
  if(is.null(subset)){subset = m$tophits}
  return(m$Bhat[subset,])
}
