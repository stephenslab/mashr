test_that("posterior calculations with common cov match regular version",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(1,1,1))
  data = set_mash_data(Bhat,Shat)
  alpha=0
  if (alpha != 0 && !all(data$Shat == 1)) {
    ## alpha models dependence of effect size on standard error
    ## alpha > 0 implies larger effects has large standard error
    ## a special case when alpha = 1 is the EZ model
    data$Shat_alpha = data$Shat^alpha
    data$Bhat = data$Bhat / data$Shat_alpha
    data$Shat = data$Shat^(1-alpha)
  } else {
    data$Shat_alpha = matrix(1, nrow(data$Shat), ncol(data$Shat))
  }
  data$alpha = alpha
  Ulist = cov_canonical(data)
  posterior_weights = matrix(1/length(Ulist),nrow = 2, ncol=length(Ulist))
  out1 = compute_posterior_matrices_general_R(data,A=diag(3),Ulist,posterior_weights)
  out2 = compute_posterior_matrices_common_cov_R(data,A=diag(3),Ulist,posterior_weights)
  expect_equal(out1,out2)
}
)
