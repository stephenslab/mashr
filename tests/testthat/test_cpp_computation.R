test_that("compare likelihood computation R vs C++ in provided data of 100 X 5", {

  load("calc_lik_matrix_data.RData")
  ## Get the number of samples (J) and the number of mixture components (P).
  J <- nrow(data$Bhat)
  P <- length(Ulist)

  ## Compute the likelihood matrix using the R implementation.
  cat(sprintf("Computing %d x %d likelihood matrix using R version.\n",J,P))
  out.time <- system.time(out1 <- calc_lik_matrix(data,Ulist,log = TRUE,
                                                  algorithm.version = "R"))
  cat(sprintf("Likelihood calculations took %0.2f seconds.\n",
              out.time["elapsed"]))
  
  ## Compute the likelihood matrix using the Rcpp implementation.
  cat(sprintf("Computing %d x %d likelihood matrix using Rcpp version.\n",J,P))
  out.time <- system.time(out2 <- calc_lik_matrix(data,Ulist,log = TRUE,
                                                  algorithm.version = "Rcpp"))
  cat(sprintf("Likelihood calculations took %0.2f seconds.\n",
              out.time["elapsed"]))
  
  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare likelihood computation R vs C++ in simulated data R = 1", {
  set.seed(999)
  N = 10
  Shat = matrix(rep(1,N),ncol=1)
  Bhat = matrix(rnorm(N),ncol=1)
  data = set_mash_data(Bhat,Shat)
  Ulist = list(id=matrix(1,nrow=1))
  out1 = calc_lik_matrix(data,Ulist)
  out2 = calc_lik_matrix(data,Ulist,algorithm.version = "R")
  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare likelihood computation R vs C++ in simulated data common cov", {
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(1,1,1))
  data = set_mash_data(Bhat, Shat)
  Ulist = cov_singletons(data)
  out1 = calc_lik_matrix(data, Ulist)
  out2 = calc_lik_matrix(data, Ulist, algorithm.version = "R")
  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare posterior computation R vs C++ in provided data of 100 X 5", {
  ## Load the data.
  load("compute_posterior_matrices_data.RData")
  ## Compute the posterior quantities using the R implementation.
  cat("Computing posterior matrices using R version.\n")
  out.time <- system.time(out1 <-
                compute_posterior_matrices(data,Ulist,posterior_weights,
                                           algorithm.version = "R"))
  cat(sprintf("Computation took %0.2f seconds.\n",out.time["elapsed"]))

  ## Compute the posterior quantities using the Rcpp implementation.
  cat("Computing posterior matrices using Rcpp version.\n")
  out.time <- system.time(out2 <-
                compute_posterior_matrices(data,Ulist,posterior_weights,
                                           algorithm.version = "Rcpp"))
  cat(sprintf("Computation took %0.2f seconds.\n",out.time["elapsed"]))
  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare posterior computation R vs C++ in simulated data R = 1", {
  set.seed(999)
  N = 10
  Shat = matrix(rep(1,N),ncol=1)
  Bhat = matrix(rnorm(N),ncol=1)
  data = set_mash_data(Bhat,Shat)
  Ulist = list(id=matrix(1,nrow=1))
  posterior_weights = matrix(rep(1,N),N,1)
  out1 <- compute_posterior_matrices(data,Ulist,posterior_weights,
                                     algorithm.version = "R")
  out2 <- compute_posterior_matrices(data,Ulist,posterior_weights,
                                     algorithm.version = "Rcpp")
  expect_equal(out1, out2, tolerance=1e-5)
}
)
