test_that("get functions work",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = set_mash_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  m = mash(data,Ulist,grid=c(0.5,1,2))
  expect_length(get_lfsr(m),6)
  expect_length(get_lfdr(m),6)
  expect_length(get_pm(m),6)
  expect_length(get_psd(m),6)
  expect_length(get_np(m),6)
}
)
