#' @export
estimate_null_correlation_mle_em = function(data, Ulist, init, max_iter = 10, tol=1e-4, grad = TRUE,
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
  m.model = mash(data, Ulist, prior=prior, verbose = FALSE, outputlevel = 3)
  xUlist = expand_cov(m.model$fitted_g$Ulist, m.model$fitted_g$grid, m.model$fitted_g$usepointmass)
  pi_s = get_estimated_pi(m.model, dimension = 'all')
  prior.v <- set_prior(length(pi_s), prior)

  # compute loglikelihood
  log_liks = rep(NA,max_iter+1)
  log_liks[1] = get_loglik(m.model)
  tracking = list()

  V = init

  result = list(V = V, mash.model = m.model)

  lump = list(data = data, model = m.model)

  for (i in 1:max_iter){
    if (track_fit)
      tracking[[i]] = result

    # max_V
    theta = cor_theta(V)
    if(grad){
      if(R!=2){
        stop('Gradient only available for R=2')
      }
      print('Using gradient')
      opt = optim(theta, neglogLik_m, gr = gradient, input = lump)
    }else{
      print('No gradient')
      opt = nlminb(theta, neglogLik_m, input = lump)
    }
    opt = nlminb(theta, neglogLik_m, gradient = gradient, input = lump)
    V = theta_cor(opt$par, R)

    # log_liks = c(log_liks, -opt$objective)

    data = mash_update_data(data, V = V)
    m.model = mash(data, Ulist, prior=prior, verbose = FALSE, outputlevel = 3)
    # pi_s = get_estimated_pi(m.model, dimension = 'all')
    result = list(V = V, mash.model = m.model)

    log_liks[i+1] = get_loglik(m.model)
    lump$model = m.model

    if(abs(log_liks[i+1]-log_liks[i])<tol) break
  }
  log_liks = log_liks[1:(i+1)]
  result$loglik = log_liks
  if (track_fit)
    result$trace = tracking

  return(result)
}

tr = function(X){
  sum(diag(X))
}

neglogLik_m = function(theta, input){
  data = input$data
  model = input$model
  R = n_conditions(data)
  J = n_effects(data)
  V = theta_cor(theta, R)
  Vinv = solve(V)
  Eb = E_V(data, model)

  0.5*log(det(V))*J + 0.5*tr(Vinv %*% Eb)
}

gradient = function(theta, input){
  data = input$data
  model = input$model
  l = pi*exp(theta)/(exp(theta)+1)
  cosl = cos(l)
  sinl = sin(l)
  cotl = cosl/sinl
  cscl = 1/sinl
  Eb = E_V(data, model)
  dump = (exp(theta) + 1)^2

  pi*cotl/(dump) - pi*cotl*cscl^2/(dump) * (Eb[1,1]+Eb[2,2]) + pi*cscl*(cotl^2 + cscl^2)/dump
}

# E_V = function(data, m.model){
#   n = n_effects(data)
#   Z = data$Bhat/data$Shat
#   Shat = data$Shat * data$Shat_alpha
#   post.m.shat = m.model$result$PosteriorMean / Shat
#   post.sec.shat = laply(1:n, function(i) (t(m.model$result$PosteriorCov[,,i]/Shat[i,])/Shat[i,]) +
#                           tcrossprod(post.m.shat[i,])) # nx2x2 array
#   temp1 = crossprod(Z)
#   temp2 = crossprod(post.m.shat, Z) + crossprod(Z, post.m.shat)
#   temp3 = unname(aaply(post.sec.shat, c(2,3), sum))
#
#   (temp1 - temp2 + temp3)
# }

