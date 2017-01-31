test_that("get same result as ash", {
  set.seed(100)
  test = simple_sims(err_sd = 0.5)
  # The simulation consists of equal numbers of four different types of effects:
  # null, equal among conditions, present only in first condition, independent across conditions
  ashres = ashr::ash(test$Bhat[,1],test$Shat[,1],mixcompdist="normal",outputlevel=4) # get ash results for first condition

  data = set_mash_data(test$Bhat,test$Shat)
  g = add_to_g(data,
                   cov_methods = list(sing1 = list(fn = "cov_first_singleton",args=NULL)),
                   grid = ashr::get_fitted_g(ashres)$sd
                   )
  g = add_to_g(data, "null", 1, g) # add null to g
  #issues: naming, especially of nulls..
  # and addition of g to intialize_g (does it scale everything or just new things?)

  lik = calc_relative_lik_matrix(data, g$Ulist)

  res=optimize_g(data,g,prior="nullbiased",optmethod="mixIP")

  post = compute_posterior_matrices(data, res$g_opt, res$posterior_weights)

  #plot(post$post_mean[,1],get_pm(ashres))
  expect_equal(post$post_mean[,1],get_pm(ashres))
  expect_equal(post$post_pos[,1],get_pp(ashres))
  expect_equal(post$post_neg[,1],get_np(ashres))
  expect_equal(post$post_zero[,1],get_lfdr(ashres))
}
)
