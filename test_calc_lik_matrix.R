library(profmem)
library(mashr)
# source("R/mash.R")
# source("R/likelihoods.R")

# Load the data.
cat("Loading data.\n")
load("calc_lik_matrix_data.RData")

# Get the number of samples (J) and the number of mixture components (P).
J <- nrow(data$Bhat)
P <- length(Ulist)

# Compute the likelihood matrix using the R implementation.
cat(sprintf("Computing %d x %d likelihood matrix using R version.\n",J,P))
out.time <- system.time(out.mem <- profmem::profmem({
  out <- calc_lik_matrix(data,Ulist,log = TRUE,algorithm.version = "R")
},threshold = 1000))
cat(sprintf(paste("Likelihood calculations allocated %0.2f MB",
                  "and took %0.2f seconds.\n"),
            sum(out.mem$bytes,na.rm = TRUE)/1024^2,
            out.time["elapsed"]))

# Compute the likelihood matrix using the Rcpp implementation.
cat(sprintf("Computing %d x %d likelihood matrix using Rcpp version.\n",J,P))
out.time <- system.time(out.mem <- profmem::profmem({
  out2 <- calc_lik_matrix(data,Ulist,log = TRUE,algorithm.version = "Rcpp")
},threshold = 1000))
cat(sprintf(paste("Likelihood calculations allocated %0.2f MB",
                  "and took %0.2f seconds.\n"),
            sum(out.mem$bytes,na.rm = TRUE)/1024^2,
            out.time["elapsed"]))

# Compare the two likelihood matrices.
cat("Largest differences between likelihood calculations:\n")
print(range(out - out2))
