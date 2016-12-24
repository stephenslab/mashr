test_that("prior covariance computations look right", {
  Bhat = rbind(c(1,2),c(2,4))
  Shat = rbind(c(1,1),c(1,1))
  data = set_mash_data(Bhat,Shat)
  Ulist = compute_covs_singletons(data)
  expect_equal(Ulist, list(singleton_1=cbind(c(1,0),c(0,0)), singleton_2=cbind(c(0,0),c(0,1))))
  Ulist = compute_covs_allones(data)
  expect_equal(Ulist, list(onematrix=cbind(c(1,1),c(1,1))))
  Ulist = compute_covs_identity(data)
  expect_equal(Ulist, list(id=cbind(c(1,0),c(0,1))))
  Ulist = compute_Ulist(data, c(compute_covs_allones), Ulist)
  expect_equal(Ulist, list(id = cbind(c(1,0),c(0,1)), onematrix=cbind(c(1,1),c(1,1))))
  Ulist = scale_Ulist(Ulist, c(1,2))
  expect_equal(Ulist, list(id.1 = cbind(c(1,0),c(0,1)), onematrix.1 = cbind(c(1,1),c(1,1)),
               id.2 = cbind(c(2,0),c(0,2)), onematrix.2 = cbind(c(2,2),c(2,2))))
  Ulist = compute_covs_singletons(data)
  Ulist = compute_Ulist(data, compute_covs_allzeros, Ulist)
  Ulist = scale_Ulist(Ulist, c(3,4))
  expect_equal(Ulist, list(singleton_1.1 = cbind(c(3,0),c(0,0)),
                           singleton_2.1 = cbind(c(0,0),c(0,3)),
                           zeromatrix.1 = cbind(c(0,0),c(0,0)),
                           singleton_1.2 = cbind(c(4,0),c(0,0)),
                           singleton_2.2 = cbind(c(0,0),c(0,4)),
                           zeromatrix.2 = cbind(c(0,0),c(0,0))) )
}
)
