test_that("Initialize MASH data properly", {
  Bhat = rbind(c(1,2),c(2,6))
  Shat = rbind(c(1,1),c(2,2))
  dat1 = set_mash_data(Bhat, Shat, df = 2)
  dat1$Z = dat1$Bhat/dat1$Shat
  dat2 = set_mash_data(dat1$Z)
  expect_equal(dat2$Shat, matrix(1,2,2))
  dat3 = set_mash_data(Bhat, pval = 2 * pt(-abs(Bhat/Shat), df = 2))
  expect_equal(dat1$Shat, dat3$Shat)
})
