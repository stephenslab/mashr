test_that("compare likelihood computation R vs C++ in provided data of 100 X 5", {

  load("calc_lik_matrix_data.RData")
  ## Get the number of samples (J) and the number of mixture components (P).
  J <- nrow(data$Bhat)
  P <- length(Ulist)

  ## Compute the likelihood matrix using the R implementation.
  cat(sprintf("Computing %d x %d likelihood matrix using R version.\n",J,P))
  out.time <- system.time(out.mem <- profmem::profmem({
    out1 <- calc_lik_matrix(data,Ulist,log = TRUE,algorithm.version = "R")
  },threshold = 1000))
  cat(sprintf(paste("Likelihood calculations allocated %0.2f MB",
                    "and took %0.2f seconds.\n"),
              sum(out.mem$bytes,na.rm = TRUE)/1024^2,
              out.time["elapsed"]))
  ## Compute the likelihood matrix using the Rcpp implementation.
  cat(sprintf("Computing %d x %d likelihood matrix using Rcpp version.\n",J,P))
  out.time <- system.time(out.mem <- profmem::profmem({
    out2 <- calc_lik_matrix(data,Ulist,log = TRUE,algorithm.version = "Rcpp")
  },threshold = 1000))
  cat(sprintf(paste("Likelihood calculations allocated %0.2f MB",
                    "and took %0.2f seconds.\n"),
              sum(out.mem$bytes,na.rm = TRUE)/1024^2,
              out.time["elapsed"]))
  expect_equal(out1, out2, tolerance=1e-8)
}
)

test_that("compare likelihood computation R vs C++ in simulated data R = 1", {
  set.seed(999)
  N = 10
  Shat = matrix(rep(1,N),ncol=1)
  Bhat = matrix(rnorm(N),ncol=1)
  data = mashr::set_mash_data(Bhat,Shat)
  Ulist = list(id=matrix(1,nrow=1))
  out1 = calc_lik_matrix(data,Ulist)
  out2 = calc_lik_matrix(data,Ulist,algorithm.version = "R")
  expect_equal(out1, out2, tolerance=1e-8)
}
)
