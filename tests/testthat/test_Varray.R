context("test_Varray.R")
test_that("Input null correlation as 3 dim array in R version, alpha=0", {
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(1,0.5,1))

  V = matrix(c(1, 0.2,0.2,0.2,1,0.2,0.2,0.2,1), 3,3)
  V3D = array(V, dim = c(3,3,2))

  data = mash_set_data(Bhat,Shat, V = V, alpha=0)
  data.3D = mash_set_data(Bhat, Shat, V=V3D, alpha=0)

  U.c = cov_canonical(data)

  m = mash(data, U.c, algorithm.version = 'R', verbose = F)

  m3D = mash(data.3D, U.c, algorithm.version = 'R', verbose = F)

  expect_equal(m, m3D)
}
)

test_that("Input null correlation as 3 dim array in R version, alpha=1", {
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(1,0.5,1))

  V = matrix(c(1, 0.2,0.2,0.2,1,0.2,0.2,0.2,1), 3,3)
  V3D = array(V, dim = c(3,3,2))

  data = mash_set_data(Bhat,Shat, V = V, alpha=1)
  data.3D = mash_set_data(Bhat, Shat, V=V3D, alpha=1)

  U.c = cov_canonical(data)

  m = mash(data, U.c, algorithm.version = 'R', verbose = F)

  m3D = mash(data.3D, U.c, algorithm.version = 'R', verbose = F)

  expect_equal(m, m3D)
}
)
