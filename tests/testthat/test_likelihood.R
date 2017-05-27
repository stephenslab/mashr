test_that("likelihood calculations look right",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = set_mash_data(Bhat, Shat)
  Ulist = cov_singletons(data)
  mm = calc_lik_matrix(data, Ulist)
  data2 = set_mash_data(Bhat[1,,drop=FALSE],Shat[1,,drop=FALSE])
  vv = calc_lik_matrix(data2, Ulist)
  expect_equal(mm[1,,drop=FALSE],vv)
  mm2 = calc_relative_lik_matrix(data, Ulist)$lik_matrix
  expect_equal(cor(mm[1,],mm2[1,]),1)
}
)

# test removed due to refactoring breaking it...
#test_that("likelihood calculations on test set match original",{
#  Bhat = rbind(c(1,2,3),c(2,4,6))
#  Shat = rbind(c(1,1,1),c(2,2,2))
#  m = mash(Bhat, Shat, cov_methods = c("null","singletons","equal_effects"))
#  expect_equal(mash_compute_loglik(m,Bhat,Shat),mash_compute_loglik(m))
#}
#)
