context("test_precomputed.R")
muffled_chol = function(x, ...)
withCallingHandlers(chol(x, ...),
                    warning = function(w) {
                      if (grepl("the matrix is either rank-deficient or indefinite", w$message))
                        invokeRestart("muffleWarning")
                    } )
test_that("Precomputed precision matrix agrees with regular computations for logliklihood", {
  skip("Cannot test it due to numerical differences between chol in armadillo and R ...")
  set.seed(1)
  simdata = simple_sims(500,5,1)
  data = mash_set_data(simdata$Bhat, simdata$Shat, alpha = 0)
  U.c = cov_canonical(data)
  grid = autoselect_grid(data,sqrt(2))
  Ulist = normalize_Ulist(U.c)
  xUlist = expand_cov(Ulist,grid,TRUE)
  loglik1 = calc_lik_rcpp(t(data$Bhat),t(data$Shat),data$V,
                             matrix(0,0,0), simplify2array(xUlist),F,T)$data
  svs = data$Shat[1,] * t(data$V * data$Shat[1,])
  sigma_rooti = list()
  for (i in 1:length(xUlist)) sigma_rooti[[i]] = backsolve(muffled_chol(svs + xUlist[[i]], pivot=T), diag(nrow(svs)))
  loglik2 = calc_lik_common_rcpp(t(data$Bhat),
                                 simplify2array(sigma_rooti),
                                 F)$data
  expect_equal(loglik1, loglik2, tol=1E-4)
})