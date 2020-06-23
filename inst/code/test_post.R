library(mashr)

# For testing, load all packages listed in the "imports" field, and
# load all the functions defined in the source files.
# to overwrite the mashr package defaults
# source("R/mash.R")
# source("R/compute_covs.R")
# source("R/data2cov.R")
# source("R/ed.R")
# source("R/get_functions.R")
# source("R/likelihoods.R")
# source("R/opt.R")
# source("R/plots.R")
# source("R/posterior.R")
# source("R/set_data.R")
# source("R/simulations.R")

# SCRIPT PARAMETERS.
J  <- 1e4  # Number of samples.
K  <- 10   # Number of dimensions/conditions.
nc <- 20    # Number of CPUs ("cores") to use.

J <- ceiling(J/4)*4

# Generate the error s.d.'s.
s <- sqrt(rchisq(J*K,df = J-1)/(J-1))

# Generate a simulated data set.
mash_data   <- simple_sims(J/4,K,s)
mash_data   <- mash_set_data(mash_data$Bhat, mash_data$Shat)
Ulist       <- cov_canonical(mash_data)
Ulist       <- mashr:::expand_cov(Ulist,1:100)
weights     <- matrix(runif(length(Ulist)*J),J,length(Ulist))
weights     <- weights/rowSums(weights)

# Compute posterior quantities.
cat(sprintf("Computing posterior matrices using %d cores.\n",nc))
out.time <-
  system.time(out <- compute_posterior_matrices(mash_data,Ulist,weights,
                                                algorithm.version = "Rcpp",
                                                mc.cores = nc))
cat(sprintf("Posterior computations took %0.2f seconds.\n",
            out.time["elapsed"]))