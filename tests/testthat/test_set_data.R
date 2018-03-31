test_that("Initialize MASH data properly", {
  Bhat = rbind(c(1,2),c(2,6))
  Shat = rbind(c(1,1),c(2,2))
  dat1 = mash_set_data(Bhat, Shat, df = 2)
  dat1$Z = dat1$Bhat/dat1$Shat
  dat2 = mash_set_data(dat1$Z)
  expect_equal(dat2$Shat, matrix(1,2,2))
  dat3 = mash_set_data(Bhat, pval = 2 * pt(-abs(Bhat/Shat), df = 2))
  expect_equal(dat1$Shat, dat3$Shat)
})

test_that("Initialize MASH CONTRAST data properly", {
  Bhat = rbind(c(1,2,3),c(2,6,8))
  Shat = rbind(c(1,1,1),c(2,2,2))
  dat1 = mash_set_data(Bhat, Shat)
  L = diag(3); L[,1] = -1
  L = L[2:3,]
  dat2 = mash_set_data_contrast(dat1, L)
  expect_equal(dat2$Shat, matrix(c(sqrt(2),sqrt(2),sqrt(8),sqrt(8)),2,2, byrow = TRUE))
  expect_equal(dat2$Shat_orig, dat1$Shat)
  expect_equal(dat2$L, L)
  expect_equal(dat2$alpha, 0)
})
