test_that("simple simulations look right", {
  test = simple_sims()

  data = set_mash_data(test$Bhat,test$Shat)

  g = initialize_g(data,c("null","id","sing","all_ones"),c(0.5,1,2))
  res=optimize_g(data,g,prior="nullbiased",optmethod="mixIP")
  print_biggest_comp(res$g_opt)

  test = compute_posterior_matrices(data, res$g_opt, res$posterior_weights)

}
)
