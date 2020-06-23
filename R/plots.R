#' @title Plot metaplot for an effect based on posterior from mash
#'
#' @param m the result of a mash fit
#'
#' @param i index of the effect to plot
#'
#' @param xlab Character string specifying x-axis label.
#'
#' @param ylab Character string specifying y-axis label.
#'
#' @param ... Additional arguments passed to \code{\link[rmeta]{metaplot}}.
#'
#' @importFrom ashr get_pm get_psd
#'
#' @importFrom rmeta metaplot
#'
#' @examples
#' simdata = simple_sims(50,5,1)
#' data = mash_set_data(simdata$Bhat, simdata$Shat)
#' m = mash(data, cov_canonical(data))
#' mash_plot_meta(m,1)
#'
#' @export
#'
mash_plot_meta = function(m,i,xlab="Effect size", ylab="Condition",...){
  metaplot(get_pm(m)[i,],get_psd(m)[i,],xlab=xlab,ylab=ylab,...)
}
