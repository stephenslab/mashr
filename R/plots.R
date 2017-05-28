#' Plot metaplot for an effect based on posterior from mash
#' @param m the result of a mash fit
#' @param i index of the effect to plot
#' @export
plot_meta = function(m,i,xlab="Effect size", ylab="Condition",...){
  rmeta::metaplot(get_pm(m)[i,],get_psd(m)[i,],xlab=xlab,ylab=ylab,...)
}
