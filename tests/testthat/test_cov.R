test_that("prior covariance computations look right", {
  Bhat = rbind(c(1,2),c(2,4))
  Ulist = compute_covs_singletons(Bhat)
  expect_equal(Ulist, list(cbind(c(1,0),c(0,0)), cbind(c(0,0),c(0,1))))
  Ulist = compute_covs_allones(Bhat)
  expect_equal(Ulist, list(cbind(c(1,1),c(1,1))))
  Ulist = compute_covs_identity(Bhat)
  expect_equal(Ulist, list(cbind(c(1,0),c(0,1))))
  Ulist = compute_Ulist(Bhat, c(compute_covs_allones), Ulist)
  expect_equal(Ulist, list(cbind(c(1,0),c(0,1)), cbind(c(1,1),c(1,1))))
}
)
