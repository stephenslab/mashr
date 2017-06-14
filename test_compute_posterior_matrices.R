library(profmem)
library(mashr)

# Load the data.
cat("Loading data.\n")
load("compute_posterior_matrices_data.RData")

# Compute the posterior quantities using the R implementation.
cat("Computing posterior matrices using R version.\n")
out.time <- system.time(out.mem <- profmem::profmem({
  out <- compute_posterior_matrices(data,xUlist,posterior_weights,
                                    algorithm.version = "R")
},threshold = 1000))
cat(sprintf("Computation allocated %0.2f MB and took %0.2f s.\n",
            sum(out.mem$bytes,na.rm = TRUE)/1024^2,out.time["elapsed"]))

# Compute the posterior quantities using the Rcpp implementation.
cat("Computing posterior matrices using Rcpp version.\n")
out.time <- system.time(out.mem <- profmem::profmem({
  out2 <- compute_posterior_matrices(data,xUlist,posterior_weights,
                                    algorithm.version = "Rcpp")
},threshold = 1000))
cat(sprintf("Computation allocated %0.2f MB and took %0.2f s.\n",
            sum(out.mem$bytes,na.rm = TRUE)/1024^2,out.time["elapsed"]))

# Compare the outputs of the R and Rcpp versions.
cat("Largest differences between likelihood calculations:\n")
d <- t(sapply(as.list(names(out)),function (x) range(out[[x]] - out2[[x]])))
rownames(d) <- names(out)
print(d)
