test_that("likelihood calculations look right"){
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  Ulist = compute_covs_singletons(Bhat)
  mm = calc_lik_matrix(Bhat, Shat, Ulist)
  v = calc_lik_matrix(Bhat[1,],Shat[1,],Ulist)
  expect_equal(mm[1,],v)
}
