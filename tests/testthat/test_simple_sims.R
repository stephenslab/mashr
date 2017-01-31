test_that("simple simulations look right", {
  set.seed(100)
  test = simple_sims()
  # The simulation consists of equal numbers of four different types of effects:
  # null, equal among conditions, present only in first condition, independent across conditions

  data = set_mash_data(test$Bhat,test$Shat)

  g = initialize_g(data,c("null","id","sing","all_ones"),c(0.5,1,2))
  res=optimize_g(data,g,prior="nullbiased",optmethod="mixIP")
  print_biggest_comp(res$g_opt)

  post = compute_posterior_matrices(data, res$g_opt, res$posterior_weights)
  expect_lt(mean(abs(test$B-post$post_mean)),mean(abs(test$B-test$Bhat)))
}
)

test_that("simple simulations look right; larger error", {
  set.seed(100)
  test = simple_sims(err_sd = 0.5)
  # The simulation consists of equal numbers of four different types of effects:
  # null, equal among conditions, present only in first condition, independent across conditions

  data = set_mash_data(test$Bhat,test$Shat)

  g = initialize_g(data,c("null","id","sing","all_ones"),c(0.5,1,2))
  res=optimize_g(data,g,prior="nullbiased",optmethod="mixIP")
  #print_biggest_comp(res$g_opt)

  post = compute_posterior_matrices(data, res$g_opt, res$posterior_weights)

  #plot(test$Bhat[,1], post$post_mean[,1], col = rep(1:4,100))

  ashres = ashr::ash(test$Bhat[,1],test$Shat[,1])
  #plot(ashr::get_pm(ashres),post$post_mean[,1], col = rep(1:4,100))

  expect_lt( #check mean error is better than no shrinkage
    mean(abs(test$B-post$post_mean)),
    mean(abs(test$B-test$Bhat)))
  expect_lt( #check mean error is better than ash for condition 1
    mean(abs(test$B[,1]-post$post_mean[,1])),
    mean(abs(test$B[,1]-ashr::get_pm(ashres))))


  # lfsr= pmin(post$post_zero[,1] + post$post_neg[,1], post$post_zero[,1] + post$post_pos[,1])
  #
  # plot(lfsr, ashr::get_lfsr(ashres))
  # plot(lfsr, test$B[,1])
  # plot(ashr::get_lfsr(ashres), test$B[,1])
  # plot(cumsum((test$B[,1] !=0)[order(lfsr)]))
  # lines(cumsum((test$B[,1] !=0)[order(ashr::get_lfsr(ashres))]),col=2)

  }
)
