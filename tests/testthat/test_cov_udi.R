test_that("udi covariance computations look right", {
  set.seed(1)
  Bhat = matrix(rnorm(100),ncol=5)
  L = matrix(rnorm(25),nrow=5,ncol=5)
  V = L %*% t(L)
  data = set_mash_data(Bhat,1,V=V)
  expect_error(cov_udi(data, c("U","U","D","D","A")))
  expect_error(cov_udi(data, c("U","U","I","I","I")))
  expect_equal(names(cov_udi(data, c("U","U","I","I","D"))), "cov_udi_UUIID")
}
)
