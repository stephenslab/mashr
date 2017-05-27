test_that("get same result as ash", {
  set.seed(100)
  sim_data = mashr2::simple_sims(nsamp = 100, err_sd = 0.5)
  # The simulation consists of equal numbers of four different types of effects:
  # null, equal among conditions, present only in first condition, independent across conditions
  ashres = ashr::ash(sim_data$Bhat[,1],sim_data$Shat[,1],mixcompdist="normal",outputlevel=4) # get ash results for first condition

  data = set_mash_data(sim_data$Bhat, sim_data$Shat)
  U  = list(first_singleton = cov_first_singleton(data))
  res = mashr2::mash(data,
                     U ,
                     grid = ashr::get_fitted_g(ashres)$sd,
                     prior = "nullbiased")


  post = get_posterior_matrices(res)

  #plot(post$post_mean[,1],get_pm(ashres))
  expect_equal(post$post_mean[,1],ashr::get_pm(ashres))
  expect_equal(post$post_sd[,1],ashr::get_psd(ashres))
  expect_equal(post$post_pos[,1],ashr::get_pp(ashres))
  expect_equal(post$post_neg[,1],ashr::get_np(ashres))
  expect_equal(post$post_zero[,1],ashr::get_lfdr(ashres))

  m2 = mash_run_1by1(data)
  post.ash = get_posterior_matrices(m2)
  expect_equal(post.ash$lfsr[,1], ashr::get_lfsr(ashres))
  expect_equal(post.ash$post_mean[,1], ashr::get_pm(ashres))
  expect_equal(post.ash$post_sd[,1], ashr::get_psd(ashres))
}
)
