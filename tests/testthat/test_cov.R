test_that("prior covariance computations look right", {
  Bhat = rbind(c(1,2),c(2,4))
  Shat = rbind(c(1,1),c(1,1))
  data = set_mash_data(Bhat,Shat)
  Ulist = cov_singletons(data)
  expect_equal(Ulist, list(singleton_1=cbind(c(1,0),c(0,0)), singleton_2=cbind(c(0,0),c(0,1))))

  Ulist = cov_allones(data)
  expect_equal(Ulist, list(allones=cbind(c(1,1),c(1,1))))

  Ulist = cov_identity(data)
  expect_equal(Ulist, list(id=cbind(c(1,0),c(0,1))))

  Ulist = compute_Ulist_byname(data, "allones", Ulist)
  expect_equal(Ulist, list(id = cbind(c(1,0),c(0,1)), allones=cbind(c(1,1),c(1,1))))

  Ulist = scale_cov(Ulist, c(1,2))
  expect_equal(Ulist, list(id.1 = cbind(c(1,0),c(0,1)), allones.1 = cbind(c(1,1),c(1,1)),
               id.2 = cbind(c(2,0),c(0,2)), allones.2 = cbind(c(2,2),c(2,2))))

  # test naming of scaled Ulist
  Ulist = cov_singletons(data)
  Ulist = compute_Ulist(data, cov_allzeros, Ulist)
  Ulist = scale_cov(Ulist, c(3,4))
  expect_equal(Ulist, list(singleton_1.1 = cbind(c(3,0),c(0,0)),
                           singleton_2.1 = cbind(c(0,0),c(0,3)),
                           allzeros.1 = cbind(c(0,0),c(0,0)),
                           singleton_1.2 = cbind(c(4,0),c(0,0)),
                           singleton_2.2 = cbind(c(0,0),c(0,4)),
                           allzeros.2 = cbind(c(0,0),c(0,0))) )
}
)
