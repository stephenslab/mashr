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
#' @param \dots additional parameters passed to \code{flashr::flash}
#' @return Returns a list of covariance matrices
#' @examples
#' # See https://stephenslab.github.io/mashr/articles/flash_mash.html
#' # for an example
#'
#' @importFrom assertthat assert_that
#' @importFrom flashr flash flash_set_data
#' @importFrom softImpute softImpute
#' @export
#'
cov_flash = function(data, factors=c("default", "nonneg"), subset=NULL, remove_singleton=FALSE, tag=NULL, output_model=NULL, ...) {
  # Only keep factors with at least two values greater than 1 / sqrt(n)
  find_nonunique_effects <- function(fl) {
    thresh <- 1/sqrt(ncol(fl$fitted_values))
    vals_above_avg <- colSums(fl$ldf$f > thresh)
    nonuniq_effects <- which(vals_above_avg > 1)
    message(paste("Removing", length(vals_above_avg) - length(nonuniq_effects), "singleton effect vectors"))
    return(fl$ldf$f[, nonuniq_effects, drop = FALSE])
  }
  nonneg <- function(Y, K = 1) {
    # this is the flashr:::udv_si function
    suppressWarnings(ret <- softImpute(Y, rank.max = K, type = "als", lambda = 0))
    pos_sum = sum(ret$v[ret$v > 0])
    neg_sum = -sum(ret$v[ret$v < 0])
    if (neg_sum > pos_sum) {
      return(list(u = -ret$u, d = ret$d, v = -ret$v))
    } else
    return(ret)
  }
  # extract a subset of data
  if(is.null(subset)){subset = 1:n_effects(data)}
  factors = match.arg(factors)
  # set default parameters
  args = list(...)
  args$data = flash_set_data(as.matrix(data$Bhat[subset,]))
  if (!exists("init_fn", args)) {
    args$init_fn = factors
    if (factors == 'default') args$init_fn = "udv_si"
    if (factors == 'nonneg') args$init_fn = nonneg
  } 
  if (!exists("greedy", args)) args$greedy = T
  if (!exists("backfit", args)) args$backfit = T
  if (factors == "nonneg") {
    args$ebnm_fn = "ebnm_ash"
    args$ebnm_param = list(l = list(mixcompdist = "normal",
                               optmethod = "mixSQP"),
                           f = list(mixcompdist = "+uniform",
                               optmethod = "mixSQP"))
  }
  f = do.call(flash, args)
  if (remove_singleton) flash_factors = find_nonunique_effects(f)
  else flash_factors = as.matrix(f$ldf$f)
  if (!is.null(output_model)) saveRDS(list(model=f, factors=flash_factors), output_model)
  if (missing(tag)) tag = factors
  U.flash = list()
  U.flash[[paste0("tFLASH_", tag)]] = t(f$fitted_values) %*% f$fitted_values / nrow(f$fitted_values)
  if (ncol(flash_factors)>0) U.flash = c(U.flash, c(cov_from_factors(t(flash_factors), paste0("FLASH_", tag))))
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
