test_that("array of posteriors looks right", {
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  Ulist = compute_covs_singletons(Bhat)
  U1 = posterior_cov(diag(3),Ulist[[1]])
  mu1 = posterior_mean(Bhat[1,],diag(3),U1)
  post_arrays = compute_posterior_arrays(Bhat, Shat, Ulist)
  expect_equal(post_arrays$post_means[1,1,],as.vector(mu1))

  U1 = posterior_cov(diag(0.5^2,3),Ulist[[2]])
  mu1 = posterior_mean(Bhat[2,],diag(0.5^2,3),U1)
  expect_equal(post_arrays$post_means[2,2,],as.vector(mu1))

}
)
