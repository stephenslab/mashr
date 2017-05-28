#' Find effects that have lfsr < thresh in at least one condition
#' @param m the mash result (from joint or 1by1 analysis)
#' @param thresh indicates the threshold below which to set signals
#' @return a vector containing the indices of the significant effects
#' @export
get_significant_results = function(m, thresh = 0.05){
  top_lfsr = apply(get_lfsr(m),1,min)
  which(top_lfsr< thresh)
}
