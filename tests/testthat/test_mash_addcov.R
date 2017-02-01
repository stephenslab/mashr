test_that("adding covariance matrices to mash object works", {
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  m = mash_init(Bhat,Shat)
  expect_error(mash_add_cov_list(m, list(r1cov(c(1,1,1)))))
  expect_error(mash_add_cov_list(m, list(test=r1cov(c(1,1))))) #wrong size
  mash_add_cov_list(m, list(test=r1cov(c(1,1,1))))
  expect_equal(get_cov(m),list(test = cbind(c(1,1,1),c(1,1,1),c(1,1,1))))

  mash_add_cov(m,"identity")
  mash_add_cov_singletons(m)
  expect_equal(list_cov(m),c("test","identity","singletons_1","singletons_2","singletons_3"))
  mash_add_cov_r1(m, c(1,2,3))
  mash_add_cov_r1(m, c(3,4,5))
  mash_add_cov_r1(m, rbind(c(3,4,5),c(6,7,8)))
  expect_equal(list_cov(m),c("test","identity","singletons_1","singletons_2","singletons_3","rank1_1","rank1_2","rank1_3","rank1_4"))

  expect_error(get_expanded_cov(m)) # need to add grid
  mash_add_grid(m, c(0.1,1,2))
  expect_equal(length(get_expanded_cov(m)),1+3*9)

  }
)
