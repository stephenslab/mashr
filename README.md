# mashr: Multivariate Adaptive Shrinkage in R

[![Travis Build Status](https://travis-ci.org/stephenslab/mashr.svg?branch=master)](https://travis-ci.org/stephenslab/mashr)
[![codecov](https://codecov.io/gh/stephenslab/mashr/branch/master/graph/badge.svg)](https://codecov.io/gh/stephenslab/mashr)

This package implements methods to estimate and test many effects in
many conditions (or many effects on many outcomes).

The methods use Empirical Bayes methods to estimate patterns of
similarity among conditions, and then exploit those patterns of
similarity among conditions to improve accuracy of effect estimates.
See [Urbut et al][mash-paper] for details of the model and methods.

Note that this R package is a refactoring of the code originally used
to generate the results for the manuscript. The original package code is
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

## Citing this work

If you find the masr package or any of the source code in this
repository useful for your work, please cite:

> Sarah Urbut, Gao Wang, Peter Carbonetto and Matthew Stephens
> (2019). [Flexible statistical methods for estimating and testing effects in genomic studies with multiple conditions.][mash-paper]
> *Nature Genetics* **51**, 187-195.

## Setup

Please follow these steps to install mashr.

1. Unlike most packages available on CRAN, mashr is not precompiled,
   and therefore to install mashr you will need to make sure that your
   R installation is properly set up to compile packages with C++
   source; in particular, the C++ compiler programs supported by your
   version of R should be installed on your computer, and R should be
   correctly configured to call these compilers when installing
   packages from source. For more information, see the
   [CRAN documentation][cran-docs].
   
2. Install the [latest release][mashr-release] of the mashr package
   using [devtools][devtools]:

    ```R
    install.packages("devtools")
    devtools::install_github("stephenslab/mashr@v0.2-11")
    ```
   
   This command should automatically install any missing dependencies
   that are available from CRAN. This command should also
   automatically retrieve and install the latest version of the ashr
   package from Github. If it does not, you can install the ashr
   package separately using devtools:

   ```R
   devtools::install_github("stephens999/ashr")
   ```
   
3. By default, the `devtools::install_github` function does not build
   the vignettes. If you would like to build the vignettes as well,
   you will need to several additional packages, including
   [flashr][flashr], that are used only in the vignettes:

   ```R
   install.packages(c("kableExtra","corrplot"))
   devtools::install_github("stephenslab/flashr")
   devtools::install_github("stephenslab/mashr@v0.2-11",build_vignettes = TRUE)
   ```

## Developer notes

+ When any changes are made to `roxygen2` markup or the C++ code in
the src directory, simply run `devtools::document()` to update
the [RcppExports.cpp](src/RcppExports.cpp), the package namespaces
(see [NAMESPACE](NAMESPACE)), and the package documentation files (in
the man directory),

+ These are the R commands to build the website (make sure you are
connected to Internet while running these commands):

   ```R
   library(pkgdown)
   pkgdown::build_site(mathjax = FALSE)
   ```

[cran-docs]: https://cran.r-project.org/manuals.html
[mash-paper]: https://doi.org/10.1038/s41588-018-0268-8
[mashr-release]: https://github.com/stephenslab/mashr/releases/tag/v0.2-11
[devtools]: https://github.com/r-lib/devtools
[flashr]: https://github.com/stephenslab/flashr
