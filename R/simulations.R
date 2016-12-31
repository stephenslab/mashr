#' Create some simple simulated data for testing purposes
#' @param nsamp number of samples
#' @param ncond number of conditions
#' @param err_sd the standard deviation of the errors
#' @details The simulation consists of equal numbers of four different types of effects: null, equal among conditions, present only in first condition, independent across conditions
#' @export
simple_sims = function(nsamp = 100, ncond = 5, err_sd=0.01){
  B.id = matrix(rnorm(nsamp*ncond),nrow=nsamp, ncol=ncond) #independent effects
  b = rnorm(nsamp)
  B.all = matrix(rep(b,ncond), nrow=nsamp, ncol=ncond) # identical effects
  B.zero = matrix(0, nrow=nsamp, ncol=ncond) #null effects
  B.one = B.zero #  effects that occur only in condition 1
  b2 = rnorm(nsamp)
  B.one[,1] = b2

  B = rbind(B.id, B.all, B.zero, B.one)

  Shat = matrix(err_sd, nrow=nrow(B), ncol=ncol(B))
  E = matrix(rnorm(length(Shat), mean=0, sd=Shat), nrow=nrow(B),ncol=ncol(B))
  Bhat = B+E
  return(list(B=B,Bhat=Bhat,Shat=Shat))
}
