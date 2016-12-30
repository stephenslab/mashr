test_that("simple simulations look right", {
  nsamp = 100
  ncond = 5
  set.seed(1)
  B.id = matrix(rnorm(nsamp*ncond),nrow=nsamp, ncol=ncond) #independent effects
  b = rnorm(nsamp)
  B.all = matrix(rep(b,ncond), nrow=nsamp, ncol=ncond) # identical effects
  B.zero = matrix(0, nrow=nsamp, ncol=ncond) #null effects
  B.one = B.zero #  effects that occur only in condition 1
  b2 = rnorm(nsamp)
  B.one[,1] = b2

  B = rbind(B.id, B.all, B.zero, B.one)

  Shat = matrix(0.01, nrow=nrow(B), ncol=ncol(B))
  E = matrix(rnorm(length(Shat), mean=0, sd=Shat), nrow=nrow(B),ncol=ncol(B))
  Bhat = B+E

  data = set_mash_data(Bhat,Shat)

  g = initialize_g(data,c("null","id","sing","all_ones"),c(0.5,1,2))
  g2=optimize_g(data,g,prior="nullbiased")


}
)
