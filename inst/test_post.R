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
J <- 1000 
K <- 10 

J <- ceiling(J / 4) * 4
s <- sqrt(rchisq(J*K,df=J-1) / (J-1))
## Simulation
mash_data <- simple_sims(J/4, K, s)
mash_data$V <- diag(K)
Ulist <- cov_canonical(mash_data)
Ulist <- expand_cov(Ulist, 1:100)
weights <- matrix(runif(length(Ulist) * J), J, length(Ulist))
weights <- weights / rowSums(weights)

cat(sprintf("Computing Posterior.\n"))
out.time <- system.time(out.mem <- profmem::profmem({
  out <- compute_posterior_matrices(mash_data,Ulist,weights,algorithm.version = "Rcpp")
},threshold = 1000))
cat(sprintf(paste("Likelihood calculations allocated %0.2f MB",
                  "and took %0.2f seconds.\n"),
            sum(out.mem$bytes,na.rm = TRUE)/1024^2,
            out.time["elapsed"]))
