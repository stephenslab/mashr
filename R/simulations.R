#' @title Create some simple simulated data for testing purposes
#'
#' @param nsamp number of samples of each type
#'
#' @param ncond number of conditions
#'
#' @param err_sd the standard deviation of the errors
#'
#' @details The simulation consists of equal numbers of four different
#' types of effects: null, equal among conditions, present only in
#' first condition, independent across conditions
#'
#' @examples
#' simple_sims(100, 5)
#'
#' @importFrom stats rnorm
#'
#' @export
#'
simple_sims = function(nsamp = 100, ncond = 5, err_sd=0.01){

  B.id = matrix(rnorm(nsamp*ncond),nrow=nsamp, ncol=ncond) #independent effects

  b = rnorm(nsamp)
  B.all = matrix(rep(b,ncond), nrow=nsamp, ncol=ncond) # identical effects

  B.zero = matrix(0, nrow=nsamp, ncol=ncond) #null effects

  B.one = B.zero #  effects that occur only in condition 1
  b2 = rnorm(nsamp)
  B.one[,1] = b2

  B = rbind(B.zero, B.id, B.one, B.all)

  Shat = matrix(err_sd, nrow=nrow(B), ncol=ncol(B))
  E = matrix(rnorm(length(Shat), mean=0, sd=Shat), nrow=nrow(B),ncol=ncol(B))
  Bhat = B+E
  row_ids = paste0("effect_", 1:nrow(B))
  col_ids = paste0("condition_", 1:ncol(B))
  rownames(B) = row_ids
  colnames(B) = col_ids
  rownames(Bhat) = row_ids
  colnames(Bhat) = col_ids
  rownames(Shat) = row_ids
  colnames(Shat) = col_ids
  return(list(B=B,Bhat=Bhat,Shat=Shat))
}

#' @title Create some more simple simulated data for testing purposes
#'
#' @param nsamp number of samples of each type
#'
#' @param err_sd the standard deviation of the errors
#'
#' @details The simulation consists of five conditions with two types
#' of effecc those present (and identical) in first two conditions and
#' those present (and identical) in last three conditions
#'
#' @examples
#' simple_sims2(100, 5)
#'
#' @importFrom stats rnorm
#'
#' @export
#'
simple_sims2 = function(nsamp = 100, err_sd=0.01){

  ncond=5
  b1 = rnorm(nsamp)
  B.1 = matrix(cbind(b1,b1,0,0,0),nrow=nsamp, ncol=ncond) #independent effects
  b2 = rnorm(nsamp)
  B.2 = matrix(cbind(0,0,b2,b2,b2),nrow=nsamp, ncol=ncond) #independent effects

  B = rbind(B.1, B.2)

  Shat = matrix(err_sd, nrow=nrow(B), ncol=ncol(B))
  E = matrix(rnorm(length(Shat), mean=0, sd=Shat), nrow=nrow(B),ncol=ncol(B))
  Bhat = B+E
  row_ids = paste0("effect_", 1:nrow(B))
  col_ids = paste0("condition_", 1:ncol(B))
  rownames(B) = row_ids
  colnames(B) = col_ids
  rownames(Bhat) = row_ids
  colnames(Bhat) = col_ids
  rownames(Shat) = row_ids
  colnames(Shat) = col_ids
  return(list(B=B,Bhat=Bhat,Shat=Shat))
}