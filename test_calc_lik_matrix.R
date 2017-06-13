library(profmem)
library(mashr)
source("R/mash.R")
source("R/likelihoods.R")

# Load the data.
cat("Loading data.\n")
load("calc_lik_matrix_data.RData")

J <- nrow(data$Bhat)
P <- length(Ulist)

cat(sprintf("Computing %d x %d likelihood matrix using R version.\n",J,P))
out.time <- system.time(out.mem <- profmem::profmem({
  out <- calc_lik_matrix(data,Ulist,log = TRUE,version = "R")
},threshold = 1000))
cat(sprintf(paste("Likelihood calculations allocated %0.2f MB",
                  "and took %0.2f seconds.\n"),
            sum(out.mem$bytes,na.rm = TRUE)/1024^2,
            out.time["elapsed"]))

