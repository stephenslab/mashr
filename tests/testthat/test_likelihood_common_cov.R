context("test_likelihood_common_cov.R")
test_that("likelihood calculations with common cov look right",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(1,1,1))
  data = mash_set_data(Bhat, Shat)
  Ulist = cov_singletons(data)
  mm = calc_lik_matrix(data, Ulist)
  mm2 = calc_lik_matrix_common_cov(data, Ulist)
  expect_equal(is_common_cov_Shat(data),TRUE)
  expect_equal(mm,mm2)
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat,Shat)
  expect_equal(is_common_cov_Shat(data),FALSE)
}
)

