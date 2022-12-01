#' @title Perform PCA on data and return list of candidate covariance
#' matrices
#'
#' @param data a mash data object
#'
#' @param npc the number of PCs to use
#'
#' @param subset indices of the subset of data to use (set to NULL for
#' all data)
#'
#' @return Returns a list of covariance matrices: the npc rank-one
#' covariance matrices based on the first npc PCs, and the rank npc
#' covariance matrix
#' @examples
#' data = mash_set_data(Bhat = cbind(c(1,2),c(3,4)), Shat = cbind(c(1,1),c(1,1)))
#' cov_pca(data,2)
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
cov_pca = function(data,npc,subset = NULL){
  assert_that(npc>1)
  assert_that(npc<=n_conditions(data))
  if(is.null(subset)){subset = 1:n_effects(data)}
  res.svd = svd(data$Bhat[subset,],nv=npc,nu=npc)
  # FIXME: we need to think of for the EE case what to use for svd input: Bhat or Bhat/Shat
  f = res.svd$v
  Ulist = cov_from_factors(t(f), "PCA")
  d = diag(res.svd$d[1:npc])
  Ulist = c(Ulist, list("tPCA"= f %*% d^2 %*% t(f)/length(subset)))
  for (i in 1:length(Ulist)) {
    rownames(Ulist[[i]]) <- colnames(data$Bhat)
    colnames(Ulist[[i]]) <- colnames(data$Bhat)
  }
  return(Ulist)
}

#' @title Perform Empirical Bayes Matrix Factorization via FLASH and return a list of
#' candidate covariance matrices
#'
#' @param data a mash data object
#'
#' @param factors "default" to use \code{flashr} default function to initialize factors, currently \code{udv_si}.
#' "nonneg" to implement a non-negative constraint on the factors
#'
#' @param subset indices of the subset of data to use (set to NULL for
#' all data)
#' @param remove_singleton whether or not factors corresponding to singleton matrices should be removed from output
#' @param tag specific a tag to name the contents in the return objects. You may want to choose different tags for
#' different paramter combinations. Default is set to \code{init_fn} parameter in \code{flashr::flash}.
#' @param output_model if specified a filename, the FLASH model will be saved to that file in RDS format.
#'
#' @param greedy_param list containing additional parameters passed to \code{flashier::flash.add.greedy}
#' @param backfit_param list containing additional parameters passed to \code{flashier::flash.backfit}
#' @return Returns a list of covariance matrices
#' @examples
#' # See https://stephenslab.github.io/mashr/articles/flash_mash.html
#' # for an example
#'
#' @export
#'
cov_flash = function(data, factors=c("default", "nonneg"), subset=NULL, remove_singleton=FALSE, tag=NULL, output_model=NULL, greedy_args=list(), backfit_args=list()) {
  if (!requireNamespace("flashier",quietly = TRUE))
    stop("cov_flash requires package flashier")
  # Only keep factors with at least two values greater than 1 / sqrt(n)
  find_nonunique_effects <- function(fl) {
    thresh <- 1/sqrt(nrow(fl$F.pm))
    vals_above_avg <- colSums(abs(fl$F.pm) > thresh)
    nonuniq_effects <- which(vals_above_avg > 1)
    message(paste("Removing", length(vals_above_avg) - length(nonuniq_effects), "singleton effect vectors"))
    return(fl$F.pm[, nonuniq_effects, drop = FALSE])
  }
  nonneg_init <- function(fl) {
    ret <- flashier::init.fn.softImpute(fl)
    pos_sum = sum(ret[[2]][ret[[2]] > 0])
    neg_sum = -sum(ret[[2]][ret[[2]] < 0])
    if (neg_sum > pos_sum) {
      return(list(-ret[[1]], -ret[[2]]))
    } else
      return(ret)
  }
  # extract a subset of data
  if(is.null(subset)){subset = 1:n_effects(data)}
  f = flashier::flash.init(as.matrix(data$Bhat[subset,]), var.type = 2)
  factors = match.arg(factors)
  # set default greedy parameters
  greedy_args$flash = f
  greedy_args$Kmax = 100
  if (!exists("init.fn", greedy_args)) {
    if (factors == 'default') greedy_args$init.fn = flashier::init.fn.softImpute
    if (factors == 'nonneg') greedy_args$init.fn = nonneg
  }
  if (factors == 'default') greedy_args$ebnm.fn = ebnm::ebnm_point_normal
  if (factors == 'nonneg') greedy_args$ebnm.fn = c(
    ebnm::ebnm_point_normal, ebnm::ebnm_point_exponential
  )
  f = do.call(flashier::flash.add.greedy, greedy_args)
  backfit_args$flash = f
  f = do.call(flashier::flash.backfit, backfit_args)
  if (remove_singleton) flash_factors = find_nonunique_effects(f)
  else flash_factors = ldf(f)$F
  if (!is.null(output_model)) saveRDS(list(model=f, factors=flash_factors), output_model)
  if (missing(tag)) tag = factors
  U.flash = list()
  U.flash[[paste0("tFLASH_", tag)]] = t(fitted(f)) %*% fitted(f) / nrow(fitted(f))
  if (ncol(flash_factors)>0) U.flash = c(U.flash, c(cov_from_factors(t(flash_factors), paste0("FLASH_", tag))))
  for (i in 1:length(U.flash)) {
    rownames(U.flash[[i]]) <- colnames(data$Bhat)
    colnames(U.flash[[i]]) <- colnames(data$Bhat)
  }
  return(U.flash)
}

#' @title Perform "extreme deconvolution" (Bovy et al) on a subset of
#' the data
#'
#' @param data a mash data object
#'
#' @param Ulist_init a named list of covariance matrices to use to
#' initialize ED; default is to use matrices from PCs
#'
#' @param subset a subset of data to be used when ED is run (set to
#' NULL for all the data)
#'
#' @param algorithm algorithm to run ED
#'
#' @param ... other arguments to be passed to ED algorith, see
#' \code{\link{extreme_deconvolution}} for algorithm 'bovy', or
#' \code{\link{teem_wrapper}} for algorithm 'teem'
#'
#' @details Runs the extreme deconvolution algorithm from Bovy et al
#' (Annals of Applied Statistics) to estimate data-driven covariance
#' matrices. It can be initialized with, for example running \code{cov_pca} with,
#' say, 5 PCs.
#' @examples
#' data = mash_set_data(Bhat = cbind(c(1,2),c(3,4)), Shat = cbind(c(1,1),c(1,1)))
#' U_pca = cov_pca(data,2)
#' U_x = apply(data$Bhat, 2, function(x) x - mean(x))
#' U_xx = t(U_x) %*% U_x / nrow(U_x)
#' cov_ed(data,c(U_pca, list(xx = U_xx)))
#'
#' @export
#'
cov_ed = function(data, Ulist_init, subset = NULL,
                  algorithm=c('bovy', 'teem'), ...) {
  algorithm = match.arg(algorithm)
  if (algorithm=='bovy') {
    Ulist_ed = bovy_wrapper(data, Ulist_init, subset, ...)$Ulist
  } else {
    Ulist_ed = teem_wrapper(data, Ulist_init, subset, ...)$U
  }
  names(Ulist_ed) = make_names("ED", if(is.null(names(Ulist_init))) 1:length(Ulist_ed) else names(Ulist_init))
  Ulist_ed
}

# For a vector x, return the rank one matrix xx'.
r1cov=function(x){x %*% t(x)}

#' produce list of rank-1 covariance matrices corresponding to rows of f
#'
#' @param f a matrix of factors (each row is a factor)
#'
#' @param name a string indicating the name to use
#'
#' @return a list of rank one matrices whose kth element is f[k,]
#' f[k,]' and named name_k
#'
#' @keywords internal
#'
cov_from_factors = function(f, name){
  Ulist = list()
  for(i in 1:nrow(f)){
    Ulist = c(Ulist,list(r1cov(f[i,])))
  }
  names(Ulist) = paste0(name,"_",(1:nrow(f)))
  return(Ulist)
}
