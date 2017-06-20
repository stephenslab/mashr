# Short script to compare the runtimes of calc_lik_matrix against the
# likelihood calculations within the "mash" function.
library(mashr)

set.seed(1)

grid_min = function(Bhat,Shat){
  min(Shat)/10
}

grid_max = function(Bhat,Shat){
  if (all(Bhat^2 <= Shat^2)) {
    8 * grid_min(Bhat,Shat)
  }  else {
    2 * sqrt(max(Bhat^2 - Shat^2))
  }
}

autoselect_grid = function(data,mult){
  gmax = grid_max(data$Bhat, data$Shat)
  gmin = grid_min(data$Bhat, data$Shat)
  if (mult == 0) {
    return(c(0, gmax/2))
  }
  else {
    npoint = ceiling(log2(gmax/gmin)/log2(mult))
    return(mult^((-npoint):0) * gmax)
  }
}

Bhat   <- matrix(rnorm(10000,0,2),nrow = 100,ncol = 100)
Shat   <- matrix(1,nrow = 100,ncol = 100)
data   <- set_mash_data(Bhat,Shat)
U      <- cov_canonical(data)
grid   <- autoselect_grid(data,sqrt(2))
xUlist <- expand_cov(U,grid,TRUE)

cat("R version:\n")
print(system.time(out1 <- calc_lik_matrix(data,xUlist,
                                          algorithm.version = "R")))
cat("Rcpp version:\n")
print(system.time(out2 <- calc_lik_matrix(data,xUlist,
                                          algorithm.version = "Rcpp")))

cat("Full mash analysis:\n")
out <- mash(data,U,algorithm.version = "Rcpp")
