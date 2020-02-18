context("test_ash.R")
test_that("get same result as ash, EE model", {
  library(ashr)
  set.seed(100)
  sim_data <- mashr::simple_sims(nsamp = 100, err_sd = runif(400*5))
  rownames(sim_data$Bhat) = colnames(sim_data$Bhat) = NULL
  rownames(sim_data$Shat) = colnames(sim_data$Shat) = NULL

  # The simulation consists of equal numbers of four different types
  # of effects: null, equal among conditions, present only in first
  # condition, independent across conditions
  ashres = ash(sim_data$Bhat[,1],sim_data$Shat[,1],
      mixcompdist="normal",outputlevel=3) # get ash results for first condition

  data = mash_set_data(sim_data$Bhat, sim_data$Shat, alpha=0)
  U  = list(first_singleton = cov_first_singleton(data))
  out <- capture.output(
    res <- mash(data,U,grid = get_fitted_g(ashres)$sd,prior = "nullbiased",
                usepointmass = FALSE))

  expect_equal(ashr::get_pm(res)[,1],ashr::get_pm(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_psd(res)[,1],ashr::get_psd(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_np(res)[,1],ashr::get_np(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_lfdr(res)[,1],ashr::get_lfdr(ashres),tolerance = 5e-5)

  m2 = mash_1by1(data, alpha=0)
  expect_equal(ashr::get_lfsr(m2)[,1], ashr::get_lfsr(ashres))
  expect_equal(ashr::get_pm(m2)[,1], ashr::get_pm(ashres))
  expect_equal(ashr::get_psd(m2)[,1], ashr::get_psd(ashres))
})

test_that("get same result as ash, EZ model", {
  library(ashr)
  set.seed(100)
  sim_data = mashr::simple_sims(nsamp = 100, err_sd = runif(400*5))
  rownames(sim_data$Bhat) = colnames(sim_data$Bhat) = NULL
  rownames(sim_data$Shat) = colnames(sim_data$Shat) = NULL
  # The simulation consists of equal numbers of four different types
  # of effects: null, equal among conditions, present only in first
  # condition, independent across conditions
  ashres = expect_warning(ash(sim_data$Bhat[,1],sim_data$Shat[,1],
      mixcompdist="normal",outputlevel=3, alpha = 1)) # get ash results for first condition

  data = mash_set_data(Bhat = sim_data$Bhat, Shat = sim_data$Shat, alpha=1)
  U  = list(first_singleton = cov_first_singleton(data))
  out <- capture.output(
    res <- mash(data,U,grid = get_fitted_g(ashres)$sd,prior = "nullbiased",
                usepointmass = FALSE))

  expect_equal(ashr::get_pm(res)[,1],ashr::get_pm(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_psd(res)[,1],ashr::get_psd(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_np(res)[,1],ashr::get_np(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_lfdr(res)[,1],ashr::get_lfdr(ashres),tolerance = 5e-5)
})

test_that("get same result as ash under transformation, EE model", {
  library(ashr)
  set.seed(100)
  sim_data = mashr::simple_sims(nsamp = 100, err_sd = runif(400*5))
  rownames(sim_data$Bhat) = colnames(sim_data$Bhat) = NULL
  rownames(sim_data$Shat) = colnames(sim_data$Shat) = NULL

  # The simulation consists of equal numbers of four different types
  # of effects: null, equal among conditions, present only in first
  # condition, independent across conditions
  ashres = ash(sim_data$Bhat[,1],sim_data$Shat[,1],
               mixcompdist="normal",outputlevel=3, alpha = 0) # get ash results for first condition

  data = mash_set_data(Bhat = sim_data$Bhat, Shat = sim_data$Shat, alpha=0)
  U  = list(first_singleton = cov_first_singleton(data))
  out <- capture.output(
    res <- mash(data,U,grid = get_fitted_g(ashres)$sd,prior = "nullbiased",
                usepointmass = FALSE, outputlevel = 1))
  A = rbind(c(1,0,0,0,0))
  res$result = mash_compute_posterior_matrices(res, data, A=A)
  # print('FIXME: Rcpp not implemented')

  expect_equal(dim(ashr::get_pm(res)), c(400,1))
  expect_equal(as.numeric(ashr::get_pm(res)),ashr::get_pm(ashres),tolerance = 1e-5)
  expect_equal(as.numeric(ashr::get_psd(res)),ashr::get_psd(ashres),tolerance = 1e-5)
  expect_equal(as.numeric(ashr::get_np(res)),ashr::get_np(ashres),tolerance = 1e-5)
  expect_equal(as.numeric(ashr::get_lfdr(res)),ashr::get_lfdr(ashres),tolerance = 5e-5)
})

test_that("get same result as ash under transformation, EZ model", {
  library(ashr)
  set.seed(100)
  sim_data = mashr::simple_sims(nsamp = 100, err_sd = runif(400*5))
  rownames(sim_data$Bhat) = colnames(sim_data$Bhat) = NULL
  rownames(sim_data$Shat) = colnames(sim_data$Shat) = NULL
  # The simulation consists of equal numbers of four different types
  # of effects: null, equal among conditions, present only in first
  # condition, independent across conditions
  ashres = expect_warning(ash(sim_data$Bhat[,1],sim_data$Shat[,1],
               mixcompdist="normal",outputlevel=3, alpha = 1)) # get ash results for first condition

  data = mash_set_data(Bhat = sim_data$Bhat, Shat = sim_data$Shat, alpha=1)
  U  = list(first_singleton = cov_first_singleton(data))
  out <- capture.output(
    res <- mash(data,U,grid = get_fitted_g(ashres)$sd,prior = "nullbiased",
                usepointmass = FALSE, outputlevel = 1))
  A = rbind(c(1,0,0,0,0))
  res$result = mash_compute_posterior_matrices(res, data, A=A)
  # print('FIXME: Rcpp not implemented')

  expect_equal(dim(ashr::get_pm(res)), c(400,1))
  expect_equal(as.numeric(ashr::get_pm(res)),ashr::get_pm(ashres),tolerance = 1e-5)
  expect_equal(as.numeric(ashr::get_psd(res)),ashr::get_psd(ashres),tolerance = 1e-5)
  expect_equal(as.numeric(ashr::get_np(res)),ashr::get_np(ashres),tolerance = 1e-5)
  expect_equal(as.numeric(ashr::get_lfdr(res)),ashr::get_lfdr(ashres),tolerance = 5e-5)
})
