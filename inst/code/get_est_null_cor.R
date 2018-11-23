set.seed(1)
simdata = simple_sims(500,5,1)

data = mash_set_data(simdata$Bhat, simdata$Shat)
U.c = cov_canonical(data)

V.est = estimate_null_correlation(data, U.c, max_iter = 3, track_fit = TRUE)

saveRDS(V.est, 'tests/testthat/estimate_null_cor.rds')

set.seed(1)
simdata = simple_sims(500,5,0.5)

data = mash_set_data(simdata$Bhat, simdata$Shat, alpha=1)
U.c = cov_canonical(data)

V.est = estimate_null_correlation(data, U.c, max_iter = 3, track_fit = TRUE)

saveRDS(V.est, 'tests/testthat/estimate_null_cor_alpha.rds')
