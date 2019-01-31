context("test_likelihood.R")

test_that("likelihood calculations look right",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat, Shat)
  Ulist = cov_singletons(data)
  mm = calc_lik_matrix(data, Ulist)
  data2 = mash_set_data(Bhat[1,,drop=FALSE],Shat[1,,drop=FALSE])
  vv = calc_lik_matrix(data2, Ulist)
  expect_equal(mm[1,,drop=FALSE],vv)
  mm2 = calc_relative_lik_matrix(data, Ulist)$loglik_matrix
  expect_equal(cor(mm[1,],exp(mm2[1,])),1)
})

test_that("likelihood calculations on test set match original",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  out <- capture.output(m <- mash(data,Ulist,grid = c(0.5,1,2)))
  expect_equal(mash_compute_loglik(m,data),m$loglik)
  expect_equal(mash_compute_loglik(m$fitted_g,data),m$loglik)
})

test_that('likelihood calculations on mash contrast set',{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat,Shat)
  # L = diag(3); L[,1] = -1; L = L[2:3,]
  L = contrast_matrix(3,1)
  data1 = mash_set_data_contrast(data, L)
  Ulist = cov_canonical(data1)

  out <- capture.output(m1 <- mash(data1,Ulist,grid = c(0.5,1,2),
                        algorithm.version = 'R'))
  expect_equal(mash_compute_loglik(m1,data1, algorithm.version='R'),m1$loglik)

  out <- capture.output(m2 <- mash(data1,Ulist,grid = c(0.5,1,2),
                                   algorithm.version = 'Rcpp'))
  expect_equal(mash_compute_loglik(m2,data1, algorithm.version='Rcpp'),m2$loglik)
  expect_equal(m1$loglik, m2$loglik)
})
