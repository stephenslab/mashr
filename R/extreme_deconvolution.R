.fixfix <- function(fix, ngauss) {
    if (is.null(fix)) {
        fix <- rep(0, ngauss)
    } else if (fix == FALSE | fix == TRUE) {
        fix <- rep(as.integer(fix), ngauss)
    } else if (length(fix) != ngauss) {
        warning("Dimension of fix* input does not match data (set all to the first entry)!")
        fix <- rep(as.integer(fix[0]), ngauss)
    } else {
        fix <- as.integer(fix)
    }
    return(fix)
}

#' @title Density estimation using Gaussian mixtures in the presence
#' of noisy, heterogeneous and incomplete data
#' 
#' @description We present a general algorithm to infer a
#' d-dimensional distribution function given a set of heterogeneous,
#' noisy observations or samples. This algorithm reconstructs the
#' error-deconvolved or 'underlying' distribution function common to
#' all samples, even when the individual samples have unique error and
#' missing-data properties. The underlying distribution is modeled as
#' a mixture of Gaussians, which is completely general. Model
#' parameters are chosen to optimize a justified, scalar objective
#' function: the logarithm of the probability of the data under the
#' error-convolved model, where the error convolution is different for
#' each data point. Optimization is performed by an Expectation
#' Maximization (EM) algorithm, extended by a regularization technique
#' and 'split-and-merge' procedure. These extensions mitigate problems
#' with singularities and local maxima, which are often encountered
#' when using the EM algorithm to estimate Gaussian density mixtures.
#' 
#' @param ydata [ndata,dy] matrix of observed quantities
#' 
#' @param ycovar [ndata,dy] / [ndata,dy,dy] / [dy,dy,ndata] matrix,
#' list or 3D array of observational error covariances (if [ndata,dy]
#' then the error correlations are assumed to vanish)
#' 
#' @param xamp [ngauss] array of initial amplitudes (*not* [1,ngauss])
#' 
#' @param xmean [ngauss,dx] matrix of initial means
#' 
#' @param xcovar [ngauss,dx,dx] list of matrices of initial covariances
#' 
#' @param projection [ndata,dy,dx] list of projection matrices
#' 
#' @param weight [ndata] array of weights to be applied to the data points
#' 
#' @param logweight (bool, default=False) if True, weight is actually
#' log(weight)
#' 
#' @param fixamp (default=None) None, True/False, or list of bools
#' 
#' @param fixmean (default=None) None, True/False, or list of bools
#' 
#' @param fixcovar (default=None) None, True/False, or list of bools
#' 
#' @param tol (double, default=1.e-6) tolerance for convergence
#' 
#' @param maxiter (long, default= 10**9) maximum number of iterations to
#' perform
#' 
#' @param w (double, default=0.) covariance regularization parameter (of the
#' conjugate prior)
#' 
#' @param logfile basename for several logfiles (_c.log has output
#' from the c-routine; _loglike.log has the log likelihood path of all
#' the accepted routes, i.e. only parts which increase the likelihood
#' are included, during splitnmerge)
#' 
#' @param splitnmerge (int, default=0) depth to go down the splitnmerge path
#' 
#' @param maxsnm (Bool, default=False) use the maximum number of split 'n'
#' merge steps, K*(K-1)*(K-2)/2
#' 
#' @param likeonly (Bool, default=False) only compute the total log
#' likelihood of the data
#' 
#' @return \item{avgloglikedata}{avgloglikedata after convergence}
#' \item{xamp}{updated xamp} \item{xmean}{updated xmean}
#' \item{xcovar}{updated xcovar}
#' 
#' @author Jo Bovy, David W. Hogg, & Sam T. Roweis
#' 
#' @references Inferring complete distribution functions from noisy,
#' heterogeneous and incomplete observations Jo Bovy, David W. Hogg, & Sam T.
#' Roweis, Submitted to AOAS (2009) [arXiv/0905.2979]
#' 
#' @examples
#' ydata <-
#' c(2.62434536, 0.38824359, 0.47182825, -0.07296862, 1.86540763,
#'   -1.30153870, 2.74481176, 0.23879310, 1.31903910, 0.75062962,
#'   2.46210794, -1.06014071, 0.67758280, 0.61594565, 2.13376944,
#'   -0.09989127, 0.82757179, 0.12214158, 1.04221375, 1.58281521,
#'   -0.10061918, 2.14472371, 1.90159072, 1.50249434, 1.90085595,
#'   0.31627214, 0.87710977, 0.06423057, 0.73211192, 1.53035547,
#'   0.30833925, 0.60324647, 0.31282730, 0.15479436, 0.32875387,
#'   0.98733540, -0.11731035, 1.23441570, 2.65980218, 1.74204416,
#'   0.80816445, 0.11237104, 0.25284171, 2.69245460, 1.05080775,
#'   0.36300435, 1.19091548, 3.10025514, 1.12015895, 1.61720311,
#'   1.30017032, 0.64775015, -0.14251820, 0.65065728, 0.79110577,
#'   1.58662319, 1.83898341, 1.93110208, 1.28558733, 1.88514116,
#'   0.24560206, 2.25286816, 1.51292982, 0.70190717, 1.48851815,
#'   0.92442829, 2.13162939, 2.51981682, 3.18557541, -0.39649633,
#'   -0.44411380, 0.49553414, 1.16003707, 1.87616892, 1.31563495,
#'   -1.02220122, 0.69379599, 1.82797464, 1.23009474, 1.76201118,
#'   0.77767186, 0.79924193, 1.18656139, 1.41005165, 1.19829972,
#'   1.11900865, 0.32933771, 1.37756379, 1.12182127, 2.12948391,
#'   2.19891788, 1.18515642, 0.62471505, 0.36126959, 1.42349435,
#'   1.07734007, 0.65614632, 1.04359686, 0.37999916, 1.69803203,
#'   0.55287144, 2.22450770, 1.40349164, 1.59357852, -0.09491185,
#'   1.16938243, 1.74055645, 0.04629940, 0.73378149, 1.03261455,
#'   -0.37311732, 1.31515939, 1.84616065, 0.14048406, 1.35054598,
#'   -0.31228341, 0.96130449, -0.61577236, 2.12141771, 1.40890054,
#'   0.97538304, 0.22483838, 2.27375593, 2.96710175, -0.85798186,
#'   2.23616403, 2.62765075, 1.33801170, -0.19926803, 1.86334532,
#'   0.81907970, 0.39607937, -0.23005814, 1.55053750, 1.79280687,
#'   0.37646927, 1.52057634, -0.14434139, 1.80186103, 1.04656730,
#'   0.81343023, 0.89825413, 1.86888616, 1.75041164, 1.52946532,
#'   1.13770121, 1.07782113, 1.61838026, 1.23249456, 1.68255141,
#'   0.68988323, -1.43483776, 2.03882460, 3.18697965, 1.44136444,
#'   0.89984477, 0.86355526, 0.88094581, 1.01740941, -0.12201873,
#'   0.48290554, 0.00297317, 1.24879916, 0.70335885, 1.49521132,
#'   0.82529684, 1.98633519, 1.21353390, 3.19069973, -0.89636092,
#'   0.35308331, 1.90148689, 3.52832571, 0.75136522, 1.04366899,
#'   0.77368576, 2.33145711, 0.71269214, 1.68006984, 0.68019840,
#'   -0.27255875, 1.31354772, 1.50318481, 2.29322588, 0.88955297,
#'   0.38263794, 1.56276110, 1.24073709, 1.28066508, 0.92688730,
#'   2.16033857, 1.36949272, 2.90465871, 2.11105670, 1.65904980,
#'   -0.62743834, 1.60231928, 1.42028220, 1.81095167, 2.04444209)
#' 
#' ydata  <- matrix(ydata,length(ydata),1)
#' N      <- dim(ydata)[1]
#' ycovar <- ydata*0 + 0.01
#' xamp   <- c(0.5,0.5)
#' xmean  <- matrix(c(0.86447943, 0.67078879, 0.322681, 0.45087394),2,2)
#' xcovar <-
#'   list(matrix(c(0.03821028, 0.04014796, 0.04108113, 0.03173839),2,2),
#'        matrix(c(0.06219194, 0.09738021, 0.04302473, 0.06778009),2,2))
#' projection <- list()
#' for (i in 1:N)
#'   projection[[i]] = matrix(c(i%%2,(i+1)%%2),1,2)
#' res <- extreme_deconvolution(ydata, ycovar, xamp, xmean, xcovar,
#'          projection=projection, logfile="ExDeconDemo")
#' 
#' @export
#' 
extreme_deconvolution <- function(ydata, ycovar, xamp, xmean, xcovar, 
                                  projection = NULL, weight = NULL, 
                                  fixamp = NULL, fixmean = NULL, fixcovar = NULL, 
                                  tol = 1e-06, maxiter = 1e+09, 
                                  w = 0, logfile = NULL, splitnmerge = 0, 
                                  maxsnm = FALSE, likeonly = FALSE, logweight = FALSE) {
    ngauss <- length(xamp)
    if (typeof(ycovar) == "list") {
        tycovar <- unlist(lapply(ycovar, t))
        diagerrors <- FALSE
    } else if (length(dim(ycovar)) == 3) {
        tycovar <- as.vector(apply(ycovar, 3, t))
        diagerrors <- FALSE
    } else {
        # a matrix
        tycovar <- as.vector(t(ycovar))
        diagerrors <- TRUE
    }
    fixamp <- .fixfix(fixamp, ngauss)
    fixmean <- .fixfix(fixmean, ngauss)
    fixcovar <- .fixfix(fixcovar, ngauss)

    if (is.null(logfile)) {
        clog <- array(0)
        clog2 <- array(0)
    } else {
        clog <- charToRaw(paste(logfile, "c.log", sep = "_"))
        clog2 <- charToRaw(paste(logfile, "loglike.log", sep = "_"))
    }

    if (maxsnm) 
        splitnmerge <- ngauss * (ngauss - 1) * (ngauss - 2)/2
    if (is.null(projection)) {
        noprojection <- TRUE
        projection <- array(0)
    } else {
        noprojection <- FALSE
        projection <- unlist(lapply(projection, t))
    }
    if (is.null(weight)) {
        noweight <- TRUE
        logweights <- array(0)
    } else if (!logweight) {
        noweight <- FALSE
        logweights <- log(weight)
    } else {
        noweight <- FALSE
        logweights <- weight
    }

    res <- extreme_deconvolution_rcpp(
        ydata, 
        tycovar,
        projection,
        logweights,
        xamp,
        xmean, 
        unlist(lapply(xcovar, t)),
        fixamp,
        fixmean, 
        fixcovar, 
        tol, 
        maxiter, 
        likeonly, 
        w, 
        clog, 
        splitnmerge,
        clog2, 
        noprojection, 
        diagerrors, 
        noweight)

    start <- 1
    end <- 0
    for (i in 1:length(xcovar)) {
        end <- end + prod(dim(xcovar[[i]]))
        xcovar[[i]] <- matrix(res$xcovar[start:end], dim(xcovar[[i]]), byrow = TRUE)
        start <- end + 1
    }
    return(list(xmean = res$xmean, xamp = res$xamp, xcovar = xcovar, avgloglikedata = res$avgloglikedata))
} 
