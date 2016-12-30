#' Create a data object for mash analysis
#' @param Bhat an N by R matrix of observed estimates
#' @param Shat an N by R matrix of corresponding standard errors
#' @return a data object for passing into mash functions
#' @export
set_mash_data = function(Bhat,Shat){
  return(list(Bhat=Bhat, Shat=Shat))
}


n_conditions = function(data){ncol(data$Bhat)}

n_effects = function(data){nrow(data$Bhat)}
