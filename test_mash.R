library(mashr)
data <- get(load("death_time_cor_nonmash.rda"))
data$betahat <- data$betahat[1:2,1:1000]
data$sebetahat <- data$sebetahat[1:2,1:1000]
mash_data <- set_mash_data(t(betahat), t(sebetahat))
U.c = cov_canonical(mash_data)
print(names(U.c))
#############  canonical mash performance ###############
cat("starting the canonical mash \n")
m.c = mash(mash_data, U.c)
save(m.c, file = "mash_circadian_canonical.rda")
