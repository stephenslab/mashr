test_that("likelihood calculations look right",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat, Shat)
  Ulist = cov_singletons(data)
  mm = calc_lik_matrix(data, Ulist)
  data2 = mash_set_data(Bhat[1,,drop=FALSE],Shat[1,,drop=FALSE])
  vv = calc_lik_matrix(data2, Ulist)
  expect_equal(mm[1,,drop=FALSE],vv)
  mm2 = calc_relative_lik_matrix(data, Ulist)$lik_matrix
  expect_equal(cor(mm[1,],mm2[1,]),1)
}
)

test_that("likelihood calculations on test set match original",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  m = mash(data,Ulist,grid = c(0.5,1,2))
  expect_equal(mash_compute_loglik(m,data),m$loglik)
  expect_equal(mash_compute_loglik(m$fitted_g,data),m$loglik)
}
)

test_that('likelihood calculations on mash contrast set',{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat,Shat)
  L = diag(3); L[,1] = -1; L = L[2:3,]
  data1 = mash_set_data_contrast(data, L)
  Ulist = cov_canonical(data1)
  # for mash contrast data, the algorithm need to be R
  m1 = mash(data1,Ulist,grid = c(0.5,1,2), algorithm.version = 'R')
  expect_equal(mash_compute_loglik(m1,data1, algorithm.version='R'),m1$loglik)
})
