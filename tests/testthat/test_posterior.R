test_that("calculations on test set with fixg match original",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  m = mash(data,Ulist,grid=c(0.5,1,2))
  expect_equal(mash(data,g=get_fitted_g(m),fixg=TRUE)$result,m$result)
}
)

test_that("diag of posterior covariance matches posterior sd",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,0.5,1),c(1,1,1))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  m = mash(data,Ulist,outputlevel=3)
}
)
