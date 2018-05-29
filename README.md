# mashr: Multivariate Adaptive Shrinkage in R

[![Travis Build Status](https://travis-ci.org/stephenslab/mashr.svg?branch=master)](https://travis-ci.org/stephenslab/mashr)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/6xpn7vfe6tslm9wn?svg=true)](https://ci.appveyor.com/project/pcarbo/mashr)
[![codecov](https://codecov.io/gh/stephenslab/mashr/branch/master/graph/badge.svg)](https://codecov.io/gh/stephenslab/mashr)

*Welcome to mashr!* This package implements methods to estimate and
test many effects in many conditions (or many effects on many
outcomes).

The methods use Empirical Bayes methods to estimate patterns of
similarity among conditions, and then exploit those patterns of
similarity among conditions to improve accuracy of effect estimates.
See [Urbut et al](http://biorxiv.org/content/early/2017/05/09/096552)
for details of the model and methods. 

Note that this R package is a refactoring of the code originally used
to create results for the paper. The original package code is
[here](http://github.com/stephenslab/mashr-paper).

## Quick Start

1. Follow the setup instructions below.

2. See the [Introductory
Vignette](https://stephenslab.github.io/mashr/articles/intro_mash.html) for an
introduction to mashr.

3. Then work through the other vignettes to learn more about mashr:
[Introduction to mash: data-driven
covariances](https://stephenslab.github.io/mashr/articles/intro_mash_dd.html)
and [Simulation with non-canonical
matrices](https://stephenslab.github.io/mashr/articles/simulate_noncanon.html).

## Setup

Please follow these steps to install the [latest version of the
mashr package](https://github.com/stephenslab/mashr/releases/tag/v0.2-6):

1. In R, install these three R packages from CRAN:

   ```R
   install.packages(c("assertthat","mvtnorm","rmeta"))
   ```

2. Optionally, install the package used for memory profiling:

   ```R
   install.packages("profmem")
   ```

3. Optionally, install MOSEK and the Rmosek package, for faster
   optimization in the `ashr` package. See the
   [ashr Github repository](https://github.com/stephens999/ashr) for
   details.

4. Install the [ExtremeDeconvolution R package](https://github.com/jobovy/extreme-deconvolution#installation). Note that you will need to link to the
   [GNU Scientific Library](https://www.gnu.org/software/gsl) to
   build this package.

5. Once you have installed all these packages, you can install and
   load the [latest
   release](https://github.com/stephenslab/mashr/releases/tag/v0.2-6)
   of the `mashr` package:

   ```R
   library(devtools)
   install_github("stephenslab/mashr@v0.2-6")
   library(mashr)
   ```

   This command should automatically retrieve and install the `ashr`
   package from Github. If it does not, install ashr separately using
   devtools:

   ```R
   library(devtools)
   install_github("stephens999/ashr")
   ```

   Alternatively, if you have cloned the repository locally, you can
   install the package by following these steps:

   ```
   R CMD build mashr
   R CMD INSTALL mashr_0.2-6.tar.gz
   ```

## Notes

+ When any changes are made to `roxygen2` markup or the C++ code in
the src directory, simply run `devtools::document()` to update
the [RcppExports.cpp](src/RcppExports.cpp), the package namespaces
(see [NAMESPACE](NAMESPACE)), and the package documentation files (in
the man directory),

+ These are the R commands to build the website (make sure you are
connected to Internet while running these commands):

```R
library(pkgdown)
build_site(mathjax = FALSE)
```

## Citation

If the data or code in this repository are useful for your research
project, please cite our preprint:

S M Urbut, G Wang, M Stephens. Flexible statistical methods for
estimating and testing effects in genomic studies with multiple
conditions. *bioRxiv* doi:10.1101/096552.


