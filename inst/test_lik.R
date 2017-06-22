library(mashr)

# install.packages(repo=NULL,pkgs="mashr_0.1-18.tar.gz",type="source")
#
# For testing, load all packages listed in the "imports" field, and
# load all the functions defined in the source files.
# to overwrite the mashr package defaults.
#
# source("../R/mash.R")
# source("../R/compute_covs.R")
# source("../R/data2cov.R")
# source("../R/ed.R")
# source("../R/get_functions.R")
source("../R/likelihoods.R")
# source("../R/opt.R")
# source("../R/plots.R")
# source("../R/posterior.R")
# source("../R/set_data.R")
# source("../R/simulations.R")

# SCRIPT PARAMETERS.
J <- 1e4
K <- 10 
J <- ceiling(J/4)*4
s <- sqrt(rchisq(J*K,df = J-1)/(J-1))

# Generate a simulated data set.
mash_data   <- simple_sims(J/4,K,s)
mash_data$V <- diag(K)
Ulist       <- cov_canonical(mash_data)
Ulist       <- expand_cov(Ulist,1:100)

# Compute the J x P likelihood matrix.
cat(sprintf("Computing %d x %d likelihood matrix.\n",J,length(Ulist)))
out.time <- system.time(out <- calc_lik_matrix(mash_data,Ulist,log = TRUE,
                                               algorithm.version = "Rcpp"))
cat(sprintf("Likelihood calculations took %0.2f seconds.\n",
            out.time["elapsed"]))
