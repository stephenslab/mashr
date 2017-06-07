library(mashr)

# SCRIPT PARAMETERS.
i <- 1:1000  # Which samples to include.
k <- 1:2     # Which dimensions to analyze.

# Load the data.
data      <- get(load("death_time_cor_nonmash.rda"))
betahat   <- data$betahat[k,i]
sebetahat <- data$sebetahat[k,i]
mash_data <- set_mash_data(t(betahat),t(sebetahat))

# Generate covariance matrices for mixture components.
U.c <- cov_canonical(mash_data)
print(names(U.c))

# Run multivariate adaptive shrinkage ("mash") analysis.
cat("starting the canonical mash\n")
out <- mash(mash_data,U.c)

# Record session info.
print(sessionInfo())
