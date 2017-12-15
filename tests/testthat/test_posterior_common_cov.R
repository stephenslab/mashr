test_that("posterior calculations with common cov match regular version",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(1,1,1))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  posterior_weights = matrix(1/length(Ulist),nrow = 2, ncol=length(Ulist))
  out1 = compute_posterior_matrices_general_R(data,A=diag(3),Ulist,posterior_weights)
  out2 = compute_posterior_matrices_common_cov_R(data,A=diag(3),Ulist,posterior_weights)
  expect_equal(out1,out2)
}
)
