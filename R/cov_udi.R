#' @title Compute a list of covariance matrices corresponding to the
#' "Unassociated", "Directly associated" and "Indirectly associated"
#' models
#'
#' @param data a mash data object, eg as created by \code{mash_set_data}
#'
#' @param model a model matrix with R columns, where R is the number
#' of conditions in the data; each row should be a vector of length R
#' with elements "U","D" and "I" indicating whether each effect is
#' Unassociated, Directly associated or Indirectly associated
#'
#' @return a named list of covariance matrices
#'
#' @details If model is specified then this returns the covariance
#' matrices for those models. The default creates all possible models.
#' For a desription of the "Unassociated", "Directly associated" and
#' "Indirectly associated" models see Stephens M (2013), A unified
#' framework for Association Analysis with Multiple Related
#' Phenotypes, PloS ONE.
#'
#' @examples
#' data = mash_set_data(Bhat = cbind(c(1,2),c(3,4)), Shat = cbind(c(1,1),c(1,1)))
#' cov_udi(data)
#' cov_udi(data,c('I','D'))
#'
#' @importFrom plyr alply
#'
#' @export
#'
cov_udi = function(data, model = udi_model_matrix(n_conditions(data))) {
  if(is.vector(model)){model = matrix(model,nrow=1)}
  res = alply(model, 1, cov_udi_single, data = data)
  attributes(res) = NULL
  names(res) = names_cov_udi(model)
  return(res)
}

#' @title Computes the covariance matrix for a single UDI model
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param data a mash data object
#'
#' @param model a vector of length R of "U","D" and "I"
#'
#' @return returns a named list of one element
#'
#' @keywords internal
#'
cov_udi_single = function(data, model) {
  R = n_conditions(data)
  V = data$V

  res = matrix(0, ncol = R, nrow = R) # holds result
  D = which(model == "D")
  U = which(model == "U")
  I = which(model == "I")
  if ((length(D) + length(U) + length(I)) != R) {
    stop("model must be vector of length R with elements U, D, I")
  }
  if (length(D) == 0) {
    stop("model must have at least one direct association")
  }
  VDD = V[D, D, drop = FALSE]
  VDU = V[D, U, drop = FALSE]
  VUU = V[U, U, drop = FALSE]

  if (length(U) > 0) {
    U0 = VDD - VDU %*% chol2inv(chol(VUU)) %*% t(VDU)
  } else {
    U0 = VDD
  }
  res[D, D] = U0

  # Here Ic means Icomplement = (U,D)
  Ic = c(U, D)
  VIIc = V[I, Ic, drop = FALSE]
  VIcIc = V[Ic, Ic, drop = FALSE]
  VIcIcinv = chol2inv(chol(VIcIc))

  BD = VIcIcinv[, (1 + length(U)):(length(U) + length(D)), drop = FALSE] # just the D columns of that matrix

  res[I, D] = VIIc %*% BD %*% U0
  res[D, I] = t(res[I, D])

  res[I, I] = VIIc %*% BD %*% U0 %*% t(BD) %*% t(VIIc)
  return(res)
}

names_cov_udi = function(model){
  apply(model,1,function(x){paste0(c("cov_udi_",x),collapse="")})
}

#' @title Create a matrix whose rows contain all possible
#' combinations of the U,D,I models that are allowed.
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param R the number of conditions
#'
#' @return a matrix that is Nmodel by R with each row containing a
#' vector of "U", "D" and "I" characters corresponding to a valide
#' model. Constraint is that there is at least one "D" in each row
#'
#' @keywords internal
#'
udi_model_matrix = function(R) {
  all = expand.grid(rep(list(c("U", "D", "I")), R))
  nD = apply(all, 1, function(x) {
    sum(x == "D")
  })
  return(all[nD > 0, ])
}
