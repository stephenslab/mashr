#' @title Create simplest simulation, cj = mu 1 data used for contrast
#' analysis
#'
#' @param nsamp number of samples of each type
#'
#' @param ncond number of conditions
#'
#' @param err_sd the standard deviation of the errors
#'
#' @details There is no true deviation exists in this case
#'
#' @examples
#' sim_contrast1(100,5)
#'
#' @export
#'
sim_contrast1 = function(nsamp = 100, ncond = 5, err_sd=sqrt(0.5)){
  # generate scalar
  Cs = rnorm(nsamp, 10)
  C = matrix(rep(Cs,ncond), nrow=nsamp, ncol=ncond)
  Shat = matrix(err_sd, nrow=nrow(C), ncol=ncol(C))
  E = matrix(rnorm(length(Shat), mean=0, sd=Shat), nrow=nrow(C),ncol=ncol(C))
  Chat = C+E
  row_ids = paste0("sample_", 1:nrow(C))
  col_ids = paste0("condition_", 1:ncol(C))
  rownames(C) = row_ids
  colnames(C) = col_ids
  rownames(Chat) = row_ids
  colnames(Chat) = col_ids
  rownames(Shat) = row_ids
  colnames(Shat) = col_ids
  return(list(C=C,Chat=Chat,Shat=Shat))
}

#' @title Create simulation with signal data used for contrast
#' analysis.
#'
#' @param nsamp Number of samples of each type.
#'
#' @param ncond Number of conditions.
#'
#' @param err_sd The standard deviation of the errors.
#'
#' @details The first condition is the reference group. The deviations
#' are the difference between the subsequent conditions with the
#' reference group.  The simulation consists of 90% null deviations,
#' 10% non-nulls.  The non-nulls consist of equal numbers of three
#' different types of deviations: equal among conditions, present only
#' in the first subsequent condition, independent across conditions.
#'
#' @examples
#' sim_contrast2(100,5)
#'
#' @export
#'
sim_contrast2 = function(nsamp = 1000, ncond = 5, err_sd=sqrt(0.5)){

  # generate scalar
  Cs = rnorm(nsamp, mean=10, sd=1)
  C = matrix(rep(Cs,ncond), nrow=nsamp, ncol=ncond)
  # generate 4 cases delta
  # 90% null
  nsamp.alt = ceiling(0.1*nsamp)
  D.zero = matrix(0, nrow=nsamp-nsamp.alt, ncol=ncond-1)
  # 10% alt
  nsamp.id = floor(nsamp.alt/3)
  nsamp.all = floor(nsamp.alt/3)
  nsamp.one = nsamp.alt - nsamp.id - nsamp.all

  D.id = matrix(rnorm(nsamp.id*(ncond-1), sd=1),nrow=nsamp.id, ncol=ncond-1) #independent deviations

  d = rnorm(nsamp.all, sd=1)
  D.all = matrix(rep(d,ncond-1), nrow=nsamp.all, ncol=ncond-1) # identical deviations

  D.one = matrix(0, nrow=nsamp.one, ncol=ncond-1) #  deviation that occur only in condition 1
  d2 = rnorm(nsamp.one, sd=1)
  D.one[,1] = d2

  D = rbind(D.zero, D.id, D.all, D.one)

  C = C + cbind(D, rep(0,nsamp))

  Shat = matrix(err_sd, nrow=nrow(C), ncol=ncol(C))
  E = matrix(rnorm(length(Shat), mean=0, sd=Shat), nrow=nrow(C),ncol=ncol(C))
  Chat = C+E
  row_ids = paste0("sample_", 1:nrow(C))
  col_ids = paste0("condition_", 1:ncol(C))
  rownames(C) = row_ids
  colnames(C) = col_ids
  rownames(Chat) = row_ids
  colnames(Chat) = col_ids
  rownames(Shat) = row_ids
  colnames(Shat) = col_ids
  return(list(C=C,Chat=Chat,Shat=Shat))
}