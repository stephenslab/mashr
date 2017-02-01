test_that("prior covariance computations look right", {
  Bhat = rbind(c(1,2),c(2,4))
  Shat = rbind(c(1,1),c(1,1))
  data = set_mash_data(Bhat,Shat)

  Ulist = cov_singletons(data)
  expect_equal(Ulist, list(cbind(c(1,0),c(0,0)), cbind(c(0,0),c(0,1))))

  Ulist = cov_equal_effects(data)
  expect_equal(Ulist, cbind(c(1,1),c(1,1)))

  Ulist = cov_identity(data)
  expect_equal(Ulist, cbind(c(1,0),c(0,1)))

  Ulist = compute_cov(data, list(identity2 = list(fn=cov_identity)))
  expect_equal(Ulist, list(identity2=cbind(c(1,0),c(0,1))))

  Ulist = compute_cov(data, "id")
  expect_equal(Ulist, list(identity=cbind(c(1,0),c(0,1))))
  Ulist = compute_cov(data, "equal_effects", Ulist)
  expect_equal(Ulist, list(identity = cbind(c(1,0),c(0,1)), equal_effects=cbind(c(1,1),c(1,1))))

  Ulist = scale_cov(Ulist, c(1,2))
  expect_equal(Ulist, list(identity.1 = cbind(c(1,0),c(0,1)), equal_effects.1 = cbind(c(1,1),c(1,1)),
               identity.2 = cbind(c(4,0),c(0,4)), equal_effects.2 = cbind(c(4,4),c(4,4))))

  # test naming of scaled Ulist
  Ulist = compute_cov(data,"sing")
  Ulist = compute_cov(data, "null", Ulist)
  Ulist = scale_cov(Ulist, c(sqrt(3),sqrt(4)))
  expect_equal(Ulist, list(singletons_1.1 = cbind(c(3,0),c(0,0)),
                           singletons_2.1 = cbind(c(0,0),c(0,3)),
                           null.1 = cbind(c(0,0),c(0,0)),
                           singletons_1.2 = cbind(c(4,0),c(0,0)),
                           singletons_2.2 = cbind(c(0,0),c(0,4)),
                           null.2 = cbind(c(0,0),c(0,0))) )

  # test matching of arguments
  Ulist = compute_cov(data, c("nul","equal"))
  expect_equal(Ulist, list(null = cbind(c(0,0),c(0,0)), equal_effects = cbind(c(1,1),c(1,1))))

}
)

test_that("normalizations look right", {
  testlist = list(cbind(c(1,2),c(3,4)), cbind(c(5,6),c(7,8)))
  expect_equal(normalize_cov(testlist[[1]]), testlist[[1]]/4)
  expect_equal(normalize_Ulist(testlist), list(testlist[[1]]/4,testlist[[2]]/8))
}
)

