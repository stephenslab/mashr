test_that("prior covariance computations look right", {
  Bhat = rbind(c(1,2),c(2,4))
  Shat = rbind(c(1,1),c(1,1))
  data = set_mash_data(Bhat,Shat)

  Ulist = cov_singletons(data)
  expect_equal(Ulist, list(cbind(c(1,0),c(0,0)), cbind(c(0,0),c(0,1))))

  Ulist = cov_all_ones(data)
  expect_equal(Ulist, cbind(c(1,1),c(1,1)))

  Ulist = cov_identity(data)
  expect_equal(Ulist, cbind(c(1,0),c(0,1)))

  Ulist = compute_cov(data, list(id2 = list(fn=cov_identity)))
  expect_equal(Ulist, list(id2=cbind(c(1,0),c(0,1))))

  Ulist = compute_cov(data, "id")
  expect_equal(Ulist, list(id=cbind(c(1,0),c(0,1))))
  Ulist = compute_cov(data, "all_ones", Ulist)
  expect_equal(Ulist, list(id = cbind(c(1,0),c(0,1)), all_ones=cbind(c(1,1),c(1,1))))

  Ulist = scale_cov(Ulist, c(1,2))
  expect_equal(Ulist, list(id.1 = cbind(c(1,0),c(0,1)), all_ones.1 = cbind(c(1,1),c(1,1)),
               id.2 = cbind(c(2,0),c(0,2)), all_ones.2 = cbind(c(2,2),c(2,2))))

  # test naming of scaled Ulist
  Ulist = compute_cov(data,"sing")
  Ulist = compute_cov(data, "all_zeros", Ulist)
  Ulist = scale_cov(Ulist, c(3,4))
  expect_equal(Ulist, list(singletons_1.1 = cbind(c(3,0),c(0,0)),
                           singletons_2.1 = cbind(c(0,0),c(0,3)),
                           all_zeros.1 = cbind(c(0,0),c(0,0)),
                           singletons_1.2 = cbind(c(4,0),c(0,0)),
                           singletons_2.2 = cbind(c(0,0),c(0,4)),
                           all_zeros.2 = cbind(c(0,0),c(0,0))) )

  # test matching of arguments
  Ulist = compute_cov(data, c("all_z","all_o"))
  expect_equal(Ulist, list(all_zeros = cbind(c(0,0),c(0,0)), all_ones = cbind(c(1,1),c(1,1))))



}
)
