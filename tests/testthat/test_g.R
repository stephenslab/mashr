test_that("basic g computations look right", {
  Bhat = rbind(c(1,2),c(2,4))
  Shat = rbind(c(1,1),c(1,1))
  data = set_mash_data(Bhat,Shat)

  g = initialize_g(data,c("null","id","sing"),c(0.5,1,2))
  expect_equal(n_comp(g),12)

  g2=optimize_g(data,g,prior="nullbiased")$g_opt
  expect_equal(n_comp(g2),12)

}
)


