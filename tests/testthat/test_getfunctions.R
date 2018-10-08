context("test_getfunctions.R")

test_that("ashr get functions work",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  out <- capture.output(m <- mash(data,Ulist,grid = c(0.5,1,2)))
  expect_length(get_lfsr(m),6)
  expect_length(get_lfdr(m),6)
  expect_length(get_pm(m),6)
  expect_length(get_psd(m),6)
  expect_length(get_np(m),6)
  expect_length(get_log10bf(m),2)
  out <- capture.output(m <- mash(data,Ulist,grid=c(0.5,1,2),
                                  usepointmass = FALSE))
  expect_null(get_log10bf(m))
  Ulist = cov_udi(data,c("I","D","U"))
  out <- capture.output(m <- mash(data,Ulist,grid = 1,normalizeU = FALSE))
  Ulist2 = c(list(null= mashr:::cov_all_zeros(data)),
             cov_udi(data,c("I","D","U")))
  temp = calc_lik_matrix(data,Ulist2,log=TRUE)
  expect_equal(as.numeric(log10(exp(temp[,2]-temp[,1]))),
               as.numeric(get_log10bf(m)),tol=1e-3)
})

test_that("get_estimated_pi works",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  out <- capture.output(m1 <- mash(data,Ulist,grid = c(0.5,1,2),
                                   usepointmass = FALSE))
  out <- capture.output(m2 <- mash(data,Ulist,grid = c(0.5,1,2),
                                   usepointmass = TRUE))

  expect_length(get_estimated_pi(m1,"all"),3*length(Ulist))
  expect_length(get_estimated_pi(m1,"cov"),length(Ulist))
  expect_length(get_estimated_pi(m1,"grid"),3)

  expect_length(get_estimated_pi(m2,"all"),3*length(Ulist)+1)
  expect_length(get_estimated_pi(m2,"cov"),length(Ulist)+1)
  expect_length(get_estimated_pi(m2,"grid"),4)
})
