# mashr: Multivariate Adaptive Shrinkage in R

[![Travis Build Status](https://travis-ci.org/stephenslab/mashr.svg?branch=master)](https://travis-ci.org/stephenslab/mashr)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/6xpn7vfe6tslm9wn?svg=true)](https://ci.appveyor.com/project/pcarbo/mashr)

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

2. See the [Introductory Vignette](docs/intro_mash.html) for an
introduction to mashr.

3. Then work through the other vignettes to learn more about mashr:
[Introduction to mash: data-driven covariances](docs/intro_mash_dd.html)
and [Simulation with non-canonical matrices](docs/simulate_noncanon.html).

## Setup

These are the necessary steps to install the most recent version of
the mashr package:

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
   load the most recent version of `mashr` available on Github:

   ```R
   library(devtools)
   install_github("stephenslab/mashr")
   library(mashr)
   ```

   This command should automatically retrieve and install
   version 2.1-19 of the `ashr` package released on Github. If it does
   not, install ashr 2.1-19 separately using devtools:

   ```R
   library(devtools)
   install_github("stephens999/ashr@v2.1-19")
   ```

   Alternatively, if you have cloned the repository locally, you can
   install the package by following these steps:

   ```
   R CMD build mashr
   R CMD INSTALL mashr_0.1-14.tar.gz
   ```

## Notes

+ When any changes are made to `roxygen2` markup or the C++ code in
the [src](src) directory, simply run `devtools::document()` to update
the [RcppExports.cpp](src/RcppExports.cpp), the package namespaces
(see [NAMESPACE](NAMESPACE)), and the package documentation files (in
the [man](man) directory),

