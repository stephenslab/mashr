context("test_simple_sims.R")

test_that("simple simulations look right", {
  set.seed(100)
  test = simple_sims()
  
  # The simulation consists of equal numbers of four different types
  # of effects: null, equal among conditions, present only in first
  # condition, independent across conditions.
  data = mash_set_data(test$Bhat, test$Shat)
  U = cov_canonical(data, c("id","sing","equal_effects"))
  out <- capture.output(
    res <- mash(data, U,grid= c(0.5,1,2), prior="nullbiased"))

  expect_lt(mean(abs(test$B-ashr::get_pm(res))),mean(abs(test$B-test$Bhat)))
}
)

test_that("simple simulations look right; larger error", {
  set.seed(100)
  test = simple_sims(err_sd = 0.5)
  
  # The simulation consists of equal numbers of four different types
  # of effects: null, equal among conditions, present only in first
  # condition, independent across conditions.
  data = mash_set_data(test$Bhat, test$Shat, alpha=0) #set alpha=0 for comparison with ash
  U = cov_canonical(data, c("id","sing","equal_effects"))
  out <- capture.output(
    res <- mash(data, U,grid= c(0.5,1,2), prior="nullbiased"))

  ashres = ashr::ash(test$Bhat[,1],test$Shat[,1])

  expect_lt( #check mean error is better than no shrinkage
    mean(abs(test$B-ashr::get_pm(res))),
    mean(abs(test$B-test$Bhat)))
  expect_lt( #check mean error is better than ash for condition 1
    mean(abs(test$B[,1]-ashr::get_pm(res)[,1])),
    mean(abs(test$B[,1]-ashr::get_pm(ashres))))


  # lfsr= pmin(post$post_zero[,1] + post$post_neg[,1], post$post_zero[,1] + post$post_pos[,1])
  #
  # plot(lfsr, ashr::get_lfsr(ashres))
  # plot(lfsr, test$B[,1])
  # plot(ashr::get_lfsr(ashres), test$B[,1])
  # plot(cumsum((test$B[,1] !=0)[order(lfsr)]))
  # lines(cumsum((test$B[,1] !=0)[order(ashr::get_lfsr(ashres))]),col=2)

  }
)
