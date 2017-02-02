test_that("get same result as ash", {
  set.seed(100)
  sim_data = mashr2::simple_sims(err_sd = 0.5)
  # The simulation consists of equal numbers of four different types of effects:
  # null, equal among conditions, present only in first condition, independent across conditions
  ashres = ashr::ash(sim_data$Bhat[,1],sim_data$Shat[,1],mixcompdist="normal",outputlevel=4) # get ash results for first condition

  res = mashr2::mash(sim_data$Bhat, sim_data$Shat,
                     cov_methods = list(sing1 = list(fn = "cov_first_singleton",args=NULL)),
                     grid = ashr::get_fitted_g(ashres)$sd,
                     prior = "nullbiased")


  post = get_posterior_matrices(res)

  #plot(post$post_mean[,1],get_pm(ashres))
  expect_equal(post$post_mean[,1],get_pm(ashres))
  expect_equal(post$post_sd[,1],get_psd(ashres))
  expect_equal(post$post_pos[,1],get_pp(ashres))
  expect_equal(post$post_neg[,1],get_np(ashres))
  expect_equal(post$post_zero[,1],get_lfdr(ashres))

}
)
