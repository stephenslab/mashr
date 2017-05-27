#' Perform PCA on data and return list of candidate covariance matrices
#' @param data a mash data object
#' @param npc the number of PCs to use
#' @param subset indices of the subset of data to use (defaults to all data)
#'
#' @return Returns a list of covariance matrices: the npc rank-one covariance matrices based on the first npc PCs,
#' and the rank npc covariance matrix
data2cov_pca = function(data,npc,subset = NULL){
  assertthat::assert_that(npc>1)
  assertthat::assert_that(npc<=n_conditions(data))
  if(is.null(subset)){subset = 1:n_effects(data)}
  res.svd = svd(data$Bhat[subset,],nv=npc,nu=npc)
  message("svd currently performed on Bhat; maybe should be Bhat/Shat?")
  f = res.svd$v
  Ulist = factors2cov(t(f), "PCA")
  d = diag(res.svd$d[1:npc])
  Ulist = c(Ulist, list("tPCA"= f %*% d^2 %*% t(f)/length(subset)))
  return(Ulist)
}

#' Perform "extreme deconvolution" (Bovy et al) on a subset of the data
#' @param data a mash data object
#' @param subset a subset of data to be used when ED is run (default of NULL means all data)
#' @param Ulist_init a named list of covariance matrices to use to initialize ED; default is to use matrices from  PCs
#' @details Runs the extreme deconvolution algorithm from Bovy et al (Annals of Applied Statistics) to estimate data-driven covariance matrices
#' The default is to initialize the EM algorithm from data2cov_pca with 5 PCs
#' @export
data2cov_ed = function(data,  subset=NULL, Ulist_init=NULL){
  if(is.null(subset)){subset = 1:n_effects(data)}
  if(is.null(Ulist_init)){Ulist_init = data2cov_pca(data,5,subset)}
  Ulist_ed = ed_wrapper(data, Ulist_init, subset)$Ulist
  names(Ulist_ed) = make_names("ED", 1:length(Ulist_ed))
  Ulist_ed
}


#' produce list of rank 1 covariance matrices corresponding to rows of f
#' @param f a matrix of factors (each row is a factor)
#' @param name a string indicating the name to use
#' @return a list of rank one matrices whose kth element is f[k,] f[k,]' and named name_k
factors2cov = function(f, name){
  Ulist = list()
  for(i in 1:nrow(f)){
    Ulist = c(Ulist,list(r1cov(f[i,])))
  }
  names(Ulist) = paste0(name,"_",(1:nrow(f)))
  return(Ulist)
}
