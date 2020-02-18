context("test_set_data.R")
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

test_that("Inf, NaN and zeros in Shat in input are detected and handled", {
  # NaN in Bhat
  Bhat = rbind(c(NaN,2),c(2,6))
  Shat = rbind(c(1,1),c(2,2))
  dat1 = expect_error(mash_set_data(Bhat, Shat, df = 2))
  # NaN in Shat
  Bhat = rbind(c(1,2),c(2,6))
  Shat = rbind(c(1,NaN),c(2,2))
  dat1 = expect_error(mash_set_data(Bhat, Shat, df = 2))
  # Inf in Bhat
  Bhat = rbind(c(Inf,2),c(2,6))
  Shat = rbind(c(1,1),c(2,2))
  dat1 = expect_error(mash_set_data(Bhat, Shat, df = 2))
  # Inf in Shat
  Bhat = rbind(c(1,2),c(2,6))
  Shat = rbind(c(1,Inf),c(2,2))
  dat1 = expect_error(mash_set_data(Bhat, Shat, df = 2))
  # Zero value in Shat
  Bhat = rbind(c(1,2),c(2,6))
  Shat = rbind(c(1,0),c(2,2))
  dat1 = expect_error(mash_set_data(Bhat, Shat, df = 2))
  # Zero value in Shat and Bhat
  Bhat = rbind(c(1,0),c(2,6))
  Shat = rbind(c(1,0),c(2,2))
  dat1 = expect_error(mash_set_data(Bhat, Shat, df = 2))
})

test_that("NA values in input data are correctly handled", {
  # NA in Bhat not in Shat
  Bhat = rbind(c(NA,2),c(2,6))
  Shat = rbind(c(1,1),c(2,2))
  dat1 = expect_error(mash_set_data(Bhat, Shat))
  # NA in both Bhat and Shat
  Bhat = rbind(c(1,NA),c(2,6))
  Shat = rbind(c(1,NA),c(2,2))
  dat1 = mash_set_data(Bhat, Shat)
  expect_equal(dat1$Bhat[1,2], 0)
  expect_equal(dat1$Shat[1,2], 1E6)
  dat1 = mash_set_data(Bhat, Shat, alpha=1)
  expect_equal(dat1$Bhat[1,2], 0)
  expect_equal(dat1$Shat[1,2], 1E6)
})

test_that("Contrast matrix generate L properly", {
  L = rbind(c(-1,1,0), c(-1,0,1))
  row.names(L) = c('2-1', '3-1')
  colnames(L) = 1:3
  attr(L, 'reference') = 1
  L.out = contrast_matrix(3, 1, 1:3)
  expect_equal(L, L.out)

  L = rbind(c(1,-1,0), c(0,-1,1))
  row.names(L) = c('1-2', '3-2')
  colnames(L) = 1:3
  attr(L, 'reference') = 2
  L.out = contrast_matrix(3, 2, 1:3)
  expect_equal(L, L.out)

  L = rbind(c(1,0,-1), c(0,1,-1))
  row.names(L) = c('1-3', '2-3')
  colnames(L) = 1:3
  attr(L, 'reference') = 3
  L.out = contrast_matrix(3, 3, 1:3)
  expect_equal(L, L.out)

  L = rbind(c(2/3, -1/3, -1/3), c(-1/3, 2/3, -1/3))
  row.names(L) = c('1-mean', '2-mean')
  colnames(L) = 1:3
  attr(L, 'reference') = 'mean'
  L.out = contrast_matrix(3, 'mean', 1:3)
  expect_equal(L, L.out)
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