test_that("get same result as ash", {
  library(ashr)
  set.seed(100)
  sim_data = mashr::simple_sims(nsamp = 100, err_sd = 0.5)
  # The simulation consists of equal numbers of four different types
  # of effects: null, equal among conditions, present only in first
  # condition, independent across conditions
  ashres = ash(sim_data$Bhat[,1],sim_data$Shat[,1],
      mixcompdist="normal",outputlevel=3) # get ash results for first condition

  data = set_mash_data(sim_data$Bhat, sim_data$Shat)
  U  = list(first_singleton = cov_first_singleton(data))
  res = mash(data,U,grid = get_fitted_g(ashres)$sd,prior = "nullbiased")

  expect_equal(ashr::get_pm(res)[,1],ashr::get_pm(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_psd(res)[,1],ashr::get_psd(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_np(res)[,1],ashr::get_np(ashres),tolerance = 1e-5)
  expect_equal(ashr::get_lfdr(res)[,1],ashr::get_lfdr(ashres),tolerance = 1e-5)

  m2 = mash_1by1(data)
  expect_equal(ashr::get_lfsr(m2)[,1], ashr::get_lfsr(ashres))
  expect_equal(ashr::get_pm(m2)[,1], ashr::get_pm(ashres))
  expect_equal(ashr::get_psd(m2)[,1], ashr::get_psd(ashres))
})
