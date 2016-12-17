compute_covs_identity = function(Bhat){
  R = ncol(Bhat)
  return(list(diag(R)))
}

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

compute_covs_allones = function(Bhat){
  R = ncol(Bhat)
  onematrix = matrix(1,nrow=R,ncol=R)
  return(list(onematrix))
}
