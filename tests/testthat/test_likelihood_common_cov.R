test_that("likelihood calculations with common cov look right",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(1,1,1))
  data = set_mash_data(Bhat, Shat)
  Ulist = cov_singletons(data)
  mm = calc_lik_matrix(data, Ulist)
  mm2 = calc_lik_matrix_common_cov(data, Ulist)
  expect_equal(mm,mm2)
}
)

test_that("posterior calculations with common cov look right",{
  set.seed(1)
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(1,1,1))
  L = matrix(rnorm(9), ncol=3)
  V = t(L) %*% L
  data = set_mash_data(Bhat, Shat, V)
  Ulist = cov_canonical(data)
  mm = calc_lik_matrix(data, Ulist)
  K = length(Ulist)
  prior = rep(1/K,K)
  post = compute_posterior_weights(prior, mm)

  expect_equal(compute_posterior_matrices(data,Ulist,post,FALSE),
               compute_posterior_matrices(data,Ulist,post,TRUE))

}
)
