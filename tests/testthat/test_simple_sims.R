test_that("simple simulations look right", {
  test = simple_sims()

  data = set_mash_data(test$Bhat,test$Shat)

  g = initialize_g(data,c("null","id","sing","all_ones"),c(0.5,1,2))
  g2=optimize_g(data,g,prior="nullbiased",optmethod="mixIP")
  print_biggest_comp(g2)

}
)
