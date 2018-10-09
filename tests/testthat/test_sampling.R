test_that("Samples from the posterior look right", {
  set.seed(100)
  test = simple_sims()
  data = mash_set_data(test$Bhat, test$Shat)
  U = cov_canonical(data)
  res = mash(data,U, algorithm.version = 'R', posterior_samples = 100)$result
  expect_equal(dim(res$PosteriorSamples), c(100,5,400))
  expect_equal(t(apply(res$PosteriorSamples, 3, colMeans)), res$PosteriorMean, tol=res$PosteriorSD)
})
