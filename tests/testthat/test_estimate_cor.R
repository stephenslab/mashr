context("test_estimate_cor.R")
test_that("Estimate null correlation properly", {
  Bhat = rbind(c(1,1.5,0.5),c(0.5,1,0.2), c(-0.8,-0.5,-1), c(0.6,1.7,0.3))
  Shat = rbind(c(1,1,1),c(1,1,1),c(1,1,1),c(1,1,1))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)

  V.adhoc = estimate_null_correlation_adhoc(data)
  data.adhoc = mash_update_data(data, V = V.adhoc)
  out <- capture.output(mash.adhoc <- mash(data.adhoc, Ulist, outputlevel = 3))

  V.est = estimate_null_correlation(data, Ulist, max_iter = 0)

  expect_equal(V.adhoc, V.est$V)
  expect_equal(mash.adhoc, V.est$mash.model)

  # saved result comes from version 0.2.18.0455
  original.null.cor = readRDS('estimate_null_cor.rds')
  set.seed(1)
  simdata = simple_sims(500,5,1)

  data = mash_set_data(simdata$Bhat, simdata$Shat)
  U.c = cov_canonical(data)

  V.est = estimate_null_correlation(data, U.c, max_iter = 3, track_fit = TRUE)
  expect_equal(V.est$V, original.null.cor$V)
  expect_equal(V.est$mash.model, original.null.cor$mash.model)
})
