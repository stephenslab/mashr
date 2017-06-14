library(profmem)
library(mashr)

# SCRIPT PARAMETERS.
i <- 1:1000  # Which samples to include. default to 1:1000
k <- 1:5     # Which dimensions to analyze. default to 1:5

# Load the data.
cat("Loading data.\n")
data      <- get(load("death_time_cor_nonmash.rda"))
betahat   <- data$betahat[k,i]
sebetahat <- data$sebetahat[k,i]
mash_data <- set_mash_data(t(betahat),t(sebetahat))

# Generate covariance matrices for mixture components.
cat("Generating covariance matrices.\n")
U.c <- cov_canonical(mash_data)
print(names(U.c))

# Run multivariate adaptive shrinkage ("mash") analysis.
cat("Running mash analysis.\n")
r <- system.time(out <- mash(mash_data,U.c,add.mem.profile = FALSE))

# Record session info.
# print(sessionInfo())
