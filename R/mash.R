# a draft of what the interface might look like


mash = function(Bhat,Shat, cov_methods = c("identity","singletons"), gridmult= sqrt(2), grid = NULL, prior="nullbiased"){
  data = set_mash_data(Bhat,Shat)
  if(is.missing(grid)){grid = autoselect_grid(data,gridmult)}

  g_init = initialize_g(data, cov_methods, grid)
  fitted_g = optimize_g(data, g_init, prior)
  result = compute_posterior_quantities(data, fitted_g)
  return(list(data=data, fitted_g= fitted_g, result=result))
}
