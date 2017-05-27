test_that("array of posteriors looks right", {
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = set_mash_data(Bhat,Shat)
  Ulist = cov_singletons(data)
  U1 = posterior_cov(diag(3),Ulist[[1]])
  mu1 = posterior_mean(Bhat[1,],diag(3),U1)
  post_array_list = compute_posterior_arrays(data, Ulist)
  expect_equal(post_array_list$post_mean[1,1,],as.vector(mu1))

  U1 = posterior_cov(diag(0.5^2,3),Ulist[[2]])
  mu1 = posterior_mean(Bhat[2,],diag(0.5^2,3),U1)
  expect_equal(post_array_list$post_mean[2,2,],as.vector(mu1))
  expect_equal(post_array_list$post_pos + post_array_list$post_neg + post_array_list$post_zero, array(1,dim=c(2,3,3)))
}
)

test_that("matrix of posterior probabilities looks right; and posterior means look right", {
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = set_mash_data(Bhat,Shat)
  Ulist = compute_cov(data,c("id","sing","equal_ef","null"))
  mm = calc_lik_matrix(data, Ulist)
  K = length(Ulist)
  prior = rep(1/K,K)
  post = compute_posterior_weights(prior, mm)
  expect_equal( cor(post[1,],mm[1,]), 1)
  expect_equal( cor(post[2,],mm[2,]), 1)

  post_array_list = compute_posterior_arrays(data, Ulist)
  post_matrix = compute_weighted_quantity(post_array_list$post_mean,post)
  expect_equal(rowSums(post*post_array_list$post_mean[,,1]),post_matrix[,1])
  expect_equal(rowSums(post*post_array_list$post_mean[,,2]),post_matrix[,2])
}
)

test_that("posterior calculations on test set match original",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = set_mash_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  m = mash(data,Ulist,grid=c(0.5,1,2))
  expect_equal(mash_compute_posterior_matrices(m,data),get_posterior_matrices(m))
}
)

