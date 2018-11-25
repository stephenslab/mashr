context("test_posterior.R")

test_that("calculations on test set with fixg match original",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,1,1),c(2,2,2))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  m <- mash(data,Ulist,grid = c(0.5,1,2), verbose = F)
  m2 <- mash(data,g = get_fitted_g(m),fixg = TRUE,verbose = F)
  expect_equal(m2$result,m$result)
})

test_that("diag of posterior covariance matches posterior sd",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,0.5,1),c(1,1,1))
  data = mash_set_data(Bhat,Shat)
  Ulist = cov_canonical(data)
  res <- mash(data,Ulist,outputlevel=3, algorithm.version = 'Rcpp', verbose = F)$result
  expect_equal(res$PosteriorSD,
               do.call(rbind,lapply(1:dim(res$PosteriorCov)[3],
                       function(i) sqrt(diag(res$PosteriorCov[,,i])))))
  res <- mash(data,Ulist,outputlevel = 3, algorithm.version = 'R', verbose = F)$result
  expect_equal(res$PosteriorSD,
               do.call(rbind,lapply(1:dim(res$PosteriorCov)[3],
                function(i) sqrt(diag(res$PosteriorCov[,,i])))))
})

test_that("diag of posterior covariance matches posterior sd: alpha=1",{
  Bhat = rbind(c(1,2,3),c(2,4,6))
  Shat = rbind(c(1,0.5,1),c(1,1,1))
  data = mash_set_data(Bhat,Shat, alpha=1)
  Ulist = cov_canonical(data)
  res <- mash(data,Ulist,outputlevel=3, algorithm.version = 'Rcpp', verbose = F)$result
  expect_equal(res$PosteriorSD,
               do.call(rbind,lapply(1:dim(res$PosteriorCov)[3],
                                    function(i) sqrt(diag(res$PosteriorCov[,,i])))))
  res <- mash(data,Ulist,outputlevel = 3, algorithm.version = 'R', verbose = F)$result
  expect_equal(res$PosteriorSD,
               do.call(rbind,lapply(1:dim(res$PosteriorCov)[3],
                                    function(i) sqrt(diag(res$PosteriorCov[,,i])))))
})
