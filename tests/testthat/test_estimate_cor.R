context("test_estimate_cor.R")
test_that("Estimate null correlation properly: alpha = 0", {
  Bhat = rbind(c(1,1.5,0.5),c(0.5,1,0.2), c(-0.8,-0.5,-1), c(0.6,1.7,0.3))
  Shat = rbind(c(1,1,1),c(1,1,1),c(1,1,1),c(1,1,1))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)

  V.adhoc = estimate_null_correlation_simple(data)
  data.adhoc = mash_update_data(data, V = V.adhoc)
  out <- capture.output(mash.adhoc <- mash(data.adhoc, Ulist, outputlevel = 3))

  V.est <- estimate_null_correlation(data, Ulist, max_iter = 0, details = TRUE)

  expect_equal(V.adhoc, V.est$V)
  expect_equal(mash.adhoc, V.est$mash.model)

  # saved result comes from version 0.2.21.0641
  skip_if_not(file.exists("estimate_null_cor.rds"))
  original.null.cor = readRDS('estimate_null_cor.rds')
  set.seed(1)
  simdata = simple_sims(500,5,0.5)

  data = mash_set_data(simdata$Bhat, simdata$Shat, alpha = 0)
  U.c = cov_canonical(data)

  V.est <- estimate_null_correlation(data, U.c, max_iter = 3, details = TRUE)
  expect_equal(V.est$V, original.null.cor$V, tolerance = 5e-4)
  expect_equal(V.est$mash.model, original.null.cor$mash.model, tolerance = 5e-4)
})

test_that("Estimate null correlation properly: alpha = 1", {
  Bhat = rbind(c(1,1.5,0.5),c(0.5,1,0.2), c(-0.8,-0.5,-1), c(0.6,1.7,0.3))
  Shat = 1.5
  data = mash_set_data(Bhat,Shat, alpha=1)
  Ulist = cov_canonical(data)

  V.adhoc = estimate_null_correlation_simple(data)
  data.adhoc = mash_update_data(data, V = V.adhoc)
  out <- capture.output(mash.adhoc <- mash(data.adhoc, Ulist, outputlevel = 3))

  V.est <- estimate_null_correlation(data, Ulist, max_iter = 0, details = TRUE)

  expect_equal(V.adhoc, V.est$V)
  expect_equal(mash.adhoc, V.est$mash.model)

  # saved result comes from version 0.2.18.0532
  skip_if_not(file.exists("estimate_null_cor_alpha.rds"))
  original.null.cor = readRDS('estimate_null_cor_alpha.rds')
  set.seed(1)
  simdata = simple_sims(500,5,0.5)

  data = mash_set_data(simdata$Bhat, simdata$Shat, alpha=1)
  U.c = cov_canonical(data)

  V.est <- estimate_null_correlation(data, U.c, max_iter = 3, details = TRUE)
  expect_equal(V.est$V, original.null.cor$V, tolerance = 5e-4)
  expect_equal(V.est$mash.model, original.null.cor$mash.model, tolerance = 5e-4)
})
