# a draft of what the interface might look like
#todo
# test posterior calculations; compare posterior means with ash estimates, especially in independent case
#  implement possible filter of data before data-driven covs?
# implement data-driven covs (pca)

mash = function(Bhat,Shat, cov_methods = c("identity","singletons"), gridmult= sqrt(2), grid = NULL, prior="nullbiased"){
  data = set_mash_data(Bhat,Shat)
  if(is.missing(grid)){grid = autoselect_grid(data,gridmult)}
  #filtered_data = filter_mash_data(data) extract top Z scores
  g_init = initialize_g(data, cov_methods, grid)
  res = optimize_g(data, g_init, prior)
  posterior_matrices = compute_posterior_matrices(data, res$g_opt, res$posterior_weights)
  return(list(data=data, fitted_g= g_opt, result=posterior_matrices))
}
