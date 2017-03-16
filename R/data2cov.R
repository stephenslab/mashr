#' Perform PCA on data and return list of candidate covariance matrices
#' @param data a mash data object
#' @param npc the number of PCs to use
#' @param subset indices of the subset of data to use (defaults to all data)
#' Returns a list of covariance matrices: the k rank-one covariance matrices based on the first k PCs,
#' and the rank k covariance matrix
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


#' produce list of rank 1 covariance matrices corresponding to rows of f
factors2cov = function(f, name){
  Ulist = list()
  for(i in 1:nrow(f)){
    Ulist = c(Ulist,list(r1cov(f[i,])))
  }
  names(Ulist) = paste0(name,"_",(1:nrow(f)))
  return(Ulist)
}
