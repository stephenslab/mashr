library(mashr)
# For testing, load all packages listed in the "imports" field, and
# load all the functions defined in the source files.
# to overwrite the mashr package defaults
library(assertthat)
library(profmem)
library(rmeta)
library(ashr)
library(mvtnorm)
source("R/mash.R")
source("R/compute_covs.R")
source("R/data2cov.R")
source("R/ed.R")
source("R/get_functions.R")
source("R/likelihoods.R")
source("R/opt.R")
source("R/plots.R")
source("R/posterior.R")
source("R/set_data.R")
source("R/simulations.R")

# SCRIPT PARAMETERS.
i <- 1:100  # Which samples to include. default to 1:1000
k <- 1:3     # Which dimensions to analyze. default to 1:5

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
r <- system.time(out <- mash(mash_data,U.c))

# Record session info.
print(sessionInfo())
