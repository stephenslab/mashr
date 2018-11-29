#' @export
estimate_null_correlation_mle = function(data, Ulist, init, max_iter = 10, tol=1e-4,
                                         track_fit=FALSE, prior = c('nullbiased', 'uniform')){
  if(class(data) != 'mash'){
    stop('data is not a "mash" object')
  }
  if(!is.null(data$L)){
    stop('We cannot estimate the null correlation for the mash contrast data.')
  }

  prior <- match.arg(prior)

  if(missing(init)){
    init = tryCatch(estimate_null_correlation_simple(data), error = function(e) FALSE)
    if(is.logical(init)){
      warning('Use Identity matrix as the initialize null correlation.')
      init = diag(n_conditions(data))
    }
  }

  R = n_conditions(data)
  data = mash_update_data(data, V = init)
  m.model = mash(data, Ulist, prior=prior, verbose = FALSE, outputlevel = 1)
  xUlist = expand_cov(m.model$fitted_g$Ulist, m.model$fitted_g$grid, m.model$fitted_g$usepointmass)
  pi_s = get_estimated_pi(m.model, dimension = 'all')
  prior.v <- set_prior(length(pi_s), prior)

  # compute loglikelihood
  log_liks = rep(NA,max_iter+1)
  log_liks[1] = get_loglik(m.model)
  tracking = list()

  V = init
  result = list(V = V, mash.model = m.model)

  lump = list(data = data, xUlist = xUlist, pi_s = pi_s)

  for (i in 1:max_iter){
    if (track_fit)
      tracking[[i]] = result

    # max_V
    theta = cor_theta(V)
    opt = nlminb(theta, neglogLik, input = lump)
    V = theta_cor(opt$par, R)

    # log_liks = c(log_liks, -opt$objective)

    data = mash_update_data(data, V = V)
    m.model = mash(data, Ulist, prior=prior, verbose = FALSE, outputlevel = 1)
    pi_s = get_estimated_pi(m.model, dimension = 'all')
    result = list(V = V, mash.model = m.model)

    log_liks[i+1] = get_loglik(m.model)
    lump$pi_s = pi_s

    if(abs(log_liks[i+1]-log_liks[i])<tol) break
  }
  log_liks = log_liks[1:(i+1)]
  result$loglik = log_liks
  if (track_fit)
    result$trace = tracking

  return(result)
}



neglogLik = function(theta, input){
  data = input$data
  xUlist = input$xUlist
  pi_s = input$pi_s

  R = n_conditions(data)
  V = theta_cor(theta, R)

  data = mash_update_data(data, V=V)

  lm = calc_relative_lik_matrix(data,xUlist)
  vloglik = compute_vloglik_from_matrix_and_pi(pi_s, lm, data$Shat_alpha)
  loglik = sum(vloglik)

  return(-loglik)
}


theta_cor = function(theta, R){

  # first get upper-triangular factor
  p = 1; q = 1
  dest = numeric((R*(R+1)/2) - 1)
  for(i in 2:R){
    aux = 1
    for(j in 2:i){
      aux1 = pi *exp(theta[p])/(1 + exp(theta[p]))
      dest[q] = aux * cos(aux1)
      aux = aux * sin(aux1)
      p = p + 1
      q = q + 1
    }
    dest[q] = aux
    q = q+1
  }

  # getting the correlations
  L = matrix(0, R,R)
  L[upper.tri(L, diag = TRUE)]= c(1,dest)

  return(crossprod(L))
}

cor_theta = function(V){
  L = chol(V)
  R = ncol(V)
  theta <- numeric(0)
  for(i in 2:R) {
    aux <- acos(L[1:(i-1),i]/sqrt(cumsum(L[i:1,i]^2)[i:2]))
    theta <- c(theta, log(aux/(pi - aux)))
  }
  return(theta)
}
