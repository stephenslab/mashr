# mashr: Multivariate Adaptive Shrinkage in R

[![Travis Build Status](https://travis-ci.org/stephenslab/mashr.svg?branch=master)](https://travis-ci.org/stephenslab/mashr)
[![codecov](https://codecov.io/gh/stephenslab/mashr/branch/master/graph/badge.svg)](https://codecov.io/gh/stephenslab/mashr)

This package implements methods to estimate and test many effects in
many conditions (or many effects on many outcomes).

The methods use Empirical Bayes methods to estimate patterns of
similarity among conditions, and then exploit those patterns of
similarity among conditions to improve accuracy of effect estimates.
See [Urbut et al][mashr-paper] for details of the model and methods.

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

## Setup

Please follow these steps to install mashr.

1. Add note here about having to compile C++ code.

2. Install the [latest release][mashr-release-latest] of the mashr
   package using [devtools][devtools]:

    ```R
    install.packages("devtools")
    library(devtools)
    devtools::install_github("stephenslab/mashr")
    ```
   
   This command should have automatically retrieved and installed the
   latest version of the ashr package from Github. If it does not,
   you can install ashr separately using devtools:

   ```R
   library(devtools)
   install_github("stephens999/ashr")
   ```

3. Optionally, install MOSEK and the Rmosek package, for faster
   optimization in the `ashr` package. See the
   [ashr Github repository](https://github.com/stephens999/ashr) for
   details.

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
build_site(mathjax = FALSE)
```

## Citation

If the data or code in this repository are useful for your research
project, please cite our preprint:

S M Urbut, G Wang, M Stephens. Flexible statistical methods for
estimating and testing effects in genomic studies with multiple
conditions. *bioRxiv* doi:10.1101/096552.

[mash-paper]: https://doi.org/10.1038/s41588-018-0268-8
[mashr-release-latest]: https://github.com/stephenslab/mashr/releases/tag/v0.2-9
[devtools]: https://github.com/r-lib/devtools
