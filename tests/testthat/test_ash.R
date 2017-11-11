test_that("get same result as ash, EE model", {
  library(ashr)
  set.seed(100)
  sim_data = mashr::simple_sims(nsamp = 100, err_sd = runif(400*5))
  rownames(sim_data$Bhat) = colnames(sim_data$Bhat) = NULL
  rownames(sim_data$Shat) = colnames(sim_data$Shat) = NULL
  # The simulation consists of equal numbers of four different types
  # of effects: null, equal among conditions, present only in first
  # condition, independent across conditions
  ashres = ash(sim_data$Bhat[,1],sim_data$Shat[,1],
      mixcompdist="normal",outputlevel=3) # get ash results for first condition

  data = set_mash_data(sim_data$Bhat, sim_data$Shat)
  U  = list(first_singleton = cov_first_singleton(data))
  res = mash(data,U,grid = get_fitted_g(ashres)$sd,prior = "nullbiased", alpha = 0, usepointmass = F)

  expect_equal(ashr::get_pm(res)[,1],ashr::get_pm(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_psd(res)[,1],ashr::get_psd(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_np(res)[,1],ashr::get_np(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_lfdr(res)[,1],ashr::get_lfdr(ashres),tolerance = 1e-5)

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
  ashres = ash(sim_data$Bhat[,1],sim_data$Shat[,1],
      mixcompdist="normal",outputlevel=3, alpha = 1) # get ash results for first condition

  data = set_mash_data(Bhat = sim_data$Bhat, Shat = sim_data$Shat)
  U  = list(first_singleton = cov_first_singleton(data))
  res = mash(data,U,grid = get_fitted_g(ashres)$sd,prior = "nullbiased", alpha=1, usepointmass = F)

  expect_equal(ashr::get_pm(res)[,1],ashr::get_pm(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_psd(res)[,1],ashr::get_psd(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_np(res)[,1],ashr::get_np(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_lfdr(res)[,1],ashr::get_lfdr(ashres),tolerance = 1e-5)
})
