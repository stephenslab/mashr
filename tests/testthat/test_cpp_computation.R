context("test_cpp_computation.R")
test_that("compare likelihood computation R vs C++ in provided data of 100 X 5", {

  load("calc_lik_matrix_data.RData")
  data$commonV = TRUE
  ## Get the number of samples (J) and the number of mixture components (P).
  J <- nrow(data$Bhat)
  P <- length(Ulist)

  ## Compute the likelihood matrix using the R implementation.
  out.time <- system.time(out1 <- calc_lik_matrix(data,Ulist,log = TRUE,
                                                  algorithm.version = "R"))

  ## Compute the likelihood matrix using the Rcpp implementation.
  out.time <- system.time(out2 <- calc_lik_matrix(data,Ulist,log = TRUE,
                                                  algorithm.version = "Rcpp",
                                                  mc.cores = 4))

  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare likelihood computation with contrast matrix R vs C++ in provided data of 100 X 5", {


  load("calc_lik_matrix_data.RData")
  data = mash_set_data(data$Bhat, data$Shat, V = data$V)
  data = mash_update_data(data, ref = 5)
  Ulist = lapply(Ulist, function(U) data$L %*% (U %*% t(data$L)))
  ## Get the number of samples (J) and the number of mixture components (P).
  J <- nrow(data$Bhat)
  P <- length(Ulist)

  ## Compute the likelihood matrix using the R implementation.
  out.time <- system.time(out1 <- calc_lik_matrix(data,Ulist,log = TRUE,
                                                  algorithm.version = "R"))

  ## Compute the likelihood matrix using the Rcpp implementation.
  out.time <- system.time(out2 <- calc_lik_matrix(data,Ulist,log = TRUE,
                                                  algorithm.version = "Rcpp",
                                                  mc.cores = 4))

  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare likelihood computation R vs C++ in simulated data R = 1", {
  set.seed(999)
  N = 10
  Shat = matrix(abs(rnorm(N)),ncol=1)
  Bhat = matrix(rnorm(N),ncol=1)
  data = mash_set_data(Bhat,Shat)
  Ulist = list(id=matrix(1,nrow=1))
  out1 = calc_lik_matrix(data,Ulist, mc.cores = 4)
  out2 = calc_lik_matrix(data,Ulist,algorithm.version = "R")
  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare likelihood computation R vs C++ in simulated data common cov", {
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(1,1,1))
  data = mash_set_data(Bhat, Shat)
  Ulist = cov_singletons(data)
  out1 = calc_lik_matrix(data, Ulist, mc.cores = 4)
  out2 = calc_lik_matrix(data, Ulist, algorithm.version = "R")
  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare likelihood computation with contrast matrix R vs C++ in simulated data common cov", {
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(1,1,1))
  data = mash_set_data(Bhat, Shat)
  data = mash_update_data(data, ref=1)
  Ulist = cov_singletons(data)
  out1 = calc_lik_matrix(data, Ulist, mc.cores = 4)
  out2 = calc_lik_matrix(data, Ulist, algorithm.version = "R")
  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare posterior computation R vs C++ in provided data of 100 X 5", {
  # Load the data.
  load("compute_posterior_matrices_data.RData")
  data$commonV = TRUE
  # Compute the posterior quantities using the R implementation.
  out.time <- system.time(out1 <-
                compute_posterior_matrices(data,Ulist,posterior_weights,
                                           algorithm.version = "R"))

  # Compute the posterior quantities using the Rcpp implementation.
  out.time <- system.time(out2 <-
                compute_posterior_matrices(data,Ulist,posterior_weights,
                                           algorithm.version = "Rcpp",
                                           mc.cores = 4))
  expect_equal(out1$PosteriorMean,out2$PosteriorMean,tolerance = 1e-5)
  expect_equal(out1$NegativeProb,out2$NegativeProb,tolerance = 1e-5)
  expect_equal(out1$lfdr,out2$lfdr,tolerance = 1e-5)
  expect_equal(out1$lfsr,out2$lfsr,tolerance = 1e-5)
}
)

test_that("compare linear transformed posterior computation R vs C++ in provided data of 100 X 5", {
  # Load the data.
  load("compute_posterior_matrices_data.RData")
  data$commonV = TRUE
  A = rbind(c(1,1,-1,-1,0),
            c(-1,1,-1,1,0))
  # Compute the posterior quantities using the R implementation.
  out.time <- system.time(out1 <-
                            compute_posterior_matrices(data,Ulist,posterior_weights,
                                                       algorithm.version = "R", A=A))

  # Compute the posterior quantities using the Rcpp implementation.
  out.time <- system.time(out2 <-
    compute_posterior_matrices(data,Ulist,posterior_weights,
                               algorithm.version = "Rcpp", A=A, mc.cores = 4))
  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare posterior computation R vs C++ in simulated data R = 1 non common cov", {
  set.seed(999)
  N = 10
  Shat = matrix(abs(rnorm(N)),ncol=1)
  Bhat = matrix(rnorm(N),ncol=1)
  data = mash_set_data(Bhat,Shat)
  Ulist = list(id=matrix(1,nrow=1))
  posterior_weights = matrix(rep(1,N),N,1)
  out1 <- compute_posterior_matrices(data,Ulist,posterior_weights,
                                     algorithm.version = "R")
  out2 <- compute_posterior_matrices(data,Ulist,posterior_weights,
                                     algorithm.version = "Rcpp",
                                     mc.cores = 4)
  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that("compare posterior computation R vs C++ in simulated data R = 1 common cov", {
  set.seed(999)
  N = 10
  Shat = matrix(rep(1,N),ncol=1)
  Bhat = matrix(rnorm(N),ncol=1)
  data = mash_set_data(Bhat,Shat)
  Ulist = list(id=matrix(1,nrow=1))
  posterior_weights = matrix(rep(1,N),N,1)
  out1 <- compute_posterior_matrices(data,Ulist,posterior_weights,
                                     algorithm.version = "R")
  out2 <- compute_posterior_matrices(data,Ulist,posterior_weights,
                                     algorithm.version = "Rcpp",
                                     mc.cores = 4)
  expect_equal(out1, out2, tolerance=1e-5)
}
)

test_that(paste("compare posterior computation R vs C++ in simulated",
                "data common cov"),{
  Bhat    <- rbind(c(1,2,3),c(2,4,6))
  Shat    <- rbind(c(1,1,1),c(1,1,1))
  data    <- mash_set_data(Bhat, Shat)
  Ulist   <- cov_canonical(data)
  Ulist   <- expand_cov(Ulist, 1:50)
  weights <- matrix(runif(length(Ulist) * 2), 2, length(Ulist))
  weights <- weights / rowSums(weights)

  out1 <- compute_posterior_matrices(data, Ulist, weights,
                                     algorithm.version = "Rcpp",
                                     mc.cores = 4)
  out2 <- compute_posterior_matrices(data, Ulist, weights,
                                     algorithm.version = "R")
  expect_equal(out1, out2, tolerance = 1e-5)
})



test_that(paste("compare posterior computation R vs C++ in simulated",
                "data non common cov"),{
  Bhat    <- rbind(c(1,2,3),c(2,4,6))
  Shat    <- rbind(c(1,2,1),c(2,1,1))
  data    <- mash_set_data(Bhat, Shat)
  Ulist   <- cov_canonical(data)
  Ulist   <- expand_cov(Ulist, 1:50)
  weights <- matrix(runif(length(Ulist) * 2), 2, length(Ulist))
  weights <- weights / rowSums(weights)

  out1 <- compute_posterior_matrices(data, Ulist, weights,
                                     algorithm.version = "Rcpp",
                                     mc.cores = 4)
  out2 <- compute_posterior_matrices(data, Ulist, weights,
                                     algorithm.version = "R")
  expect_equal(out1, out2, tolerance = 1e-5)
})

test_that(paste("compare transformed posterior computation R vs C++ in simulated",
                "data common cov"),{
  Bhat    <- rbind(c(1,2,3),c(2,4,6))
  Shat    <- rbind(c(1,1,1),c(1,1,1))
  data    <- mash_set_data(Bhat, Shat)
  Ulist   <- cov_canonical(data)
  Ulist   <- expand_cov(Ulist, 1:50)
  weights <- matrix(runif(length(Ulist) * 2), 2, length(Ulist))
  weights <- weights / rowSums(weights)

  A = rbind(c(1,1,1))

  out1 <- compute_posterior_matrices(data, Ulist, weights,
                                     algorithm.version = "Rcpp", A=A, mc.cores = 4)
  out2 <- compute_posterior_matrices(data, Ulist, weights,
                                     algorithm.version = "R", A=A)
  expect_equal(out1, out2, tolerance = 1e-5)
})

test_that(paste("compare transformed posterior computation R vs C++ in simulated",
                "data non common cov"),{
  Bhat    <- rbind(c(1,2,3),c(2,4,6))
  Shat    <- rbind(c(1,2,1),c(2,1,1))
  data    <- mash_set_data(Bhat, Shat)
  Ulist   <- cov_canonical(data)
  Ulist   <- expand_cov(Ulist, 1:50)
  weights <- matrix(runif(length(Ulist) * 2), 2, length(Ulist))
  weights <- weights / rowSums(weights)

  A = rbind(c(1,1,1))
  out1 <- compute_posterior_matrices(data, Ulist, weights,
                                     algorithm.version = "Rcpp", A=A, mc.cores = 4)
  out2 <- compute_posterior_matrices(data, Ulist, weights,
                                     algorithm.version = "R", A=A)
  expect_equal(out1, out2, tolerance = 1e-5)
})

test_that(paste("compare commonbaseline and transformed computation R vs C++ in simulated",
                "data common cov"),{
  Bhat    <- rbind(c(1,2,3),c(2,4,6))
  Shat    <- rbind(c(1,1,1),c(1,1,1))
  data    <- mash_set_data(Bhat, Shat)
  data.L  <- mash_update_data(data, ref=1)
  Ulist   <- cov_canonical(data.L)

  A = rbind(c(1,1))
  out1 <- mash(data.L, Ulist, algorithm.version = "Rcpp", A = A, verbose = F)
  out2 <- mash(data.L, Ulist, algorithm.version = "R", A = A, verbose = F)
  expect_equal(out1, out2, tolerance = 1e-5)
})

test_that(paste("compare commonbaseline and transformed computation R vs C++ in simulated",
                "data non common cov"),{
  Bhat    <- rbind(c(1,2,3),c(2,4,6))
  Shat    <- rbind(c(1,2,1),c(2,1,1))
  data    <- mash_set_data(Bhat, Shat)
  data.L  <- mash_update_data(data, ref=1)
  Ulist   <- cov_canonical(data.L)

  A = rbind(c(1,1))
  out1 <- mash(data.L, Ulist, algorithm.version = "Rcpp", A = A, verbose = F)
  out2 <- mash(data.L, Ulist, algorithm.version = "R", A = A, verbose = F)
  expect_equal(out1, out2, tolerance = 1e-5)
})

test_that(paste("compare nobaseline computation R vs C++ in simulated",
                "data common cov"),{
  Bhat    <- rbind(c(1,2,3),c(2,4,6))
  Shat    <- rbind(c(1,1,1),c(1,1,1))
  data    <- mash_set_data(Bhat, Shat)
  data.L  <- mash_update_data(data, ref='mean')
  Ulist   <- cov_canonical(data.L)

  out1 <- mash(data.L, Ulist, algorithm.version = "Rcpp", verbose = F)
  out2 <- mash(data.L, Ulist, algorithm.version = "R", verbose = F)
  expect_equal(out1, out2, tolerance = 1e-5)
})

test_that(paste("compare nobaseline computation R vs C++ in simulated",
                "data non common cov"),{
  Bhat    <- rbind(c(1,2,3),c(2,4,6))
  Shat    <- rbind(c(1,2,1),c(2,1,1))
  data    <- mash_set_data(Bhat, Shat)
  data.L  <- mash_update_data(data, ref='mean')
  Ulist   <- cov_canonical(data.L)

  out1 <- mash(data.L, Ulist, algorithm.version = "Rcpp", verbose = F)
  out2 <- mash(data.L, Ulist, algorithm.version = "R", verbose = F)
  expect_equal(out1, out2, tolerance = 1e-5)
})

test_that("Interface agree when R = 1 non common cov", {
  N = 20
  r = 1
  Shat = matrix(abs(rnorm(N)),ncol=r)
  Bhat = matrix(rnorm(N),ncol=r)
  data = mash_set_data(Bhat,Shat)

  prior = list(Ulist = list(id=diag(r)),
               grid = 1, pi_s = c(0.05,0.95), usepointmass = TRUE)
  out1 <- mash(data, g = prior, fixg = TRUE, algorithm.version = "R", verbose = F)
  out2 <- mash(data, g = prior, fixg = TRUE, algorithm.version = "Rcpp", verbose = F)
  expect_equal(get_pm(out1), get_pm(out2), tolerance=1e-5)
})

test_that("Interface agree when R > 1 non common cov, EZ version", {
  N = 20
  r = 2
  Shat = matrix(abs(rnorm(N)),ncol=r)
  Bhat = matrix(rnorm(N),ncol=r)
  data = mash_set_data(Bhat,Shat,alpha=1)

  prior = list(Ulist = list(id=diag(r)),
               grid = 1, pi_s = c(0.05,0.95), usepointmass = TRUE)
  out1 <- mash(data, g = prior, fixg = TRUE, algorithm.version = "R", verbose = F)
  out2 <- mash(data, g = prior, fixg = TRUE, algorithm.version = "Rcpp", verbose = F)
  expect_equal(get_pm(out1), get_pm(out2), tolerance=1e-5)
})
