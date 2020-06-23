# mashr: Multivariate Adaptive Shrinkage in R

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/mashr)](https://cran.r-project.org/package=mashr)
[![Travis Build Status](https://travis-ci.org/stephenslab/mashr.svg?branch=master)](https://travis-ci.org/stephenslab/mashr)
[![codecov](https://codecov.io/gh/stephenslab/mashr/branch/master/graph/badge.svg)](https://codecov.io/gh/stephenslab/mashr)

This package implements methods to estimate and test many effects in
many conditions (or many effects on many outcomes).

Empirical Bayes methods are used to estimate patterns of similarity
among conditions, then exploit these patterns of similarity among
conditions to improve accuracy of effect estimates. See
[Urbut et al][mash-paper] for details.

**Note:** the R package used to generate the results for the manuscript
is [here][mashr-pkg-for-paper].

## Quick Start

1. Install mashr from CRAN: `install.packages("mashr")`. For more
detailed setup instructions, see below.

2. See the [introductory vignette][vignette-intro] for an
introduction to mashr.

3. Then walk through these vignettes to learn more about mashr:
[Introduction to mash: data-driven covariances][vignette-data-driven-cov]
and [Simulation with non-canonical matrices][vignette-non-canonical].

## Citing this work

If you find the mashr package or any of the source code in this
repository useful for your work, please cite:

> Sarah Urbut, Gao Wang, Peter Carbonetto and Matthew Stephens
> (2019). [Flexible statistical methods for estimating and testing effects in genomic studies with multiple conditions.][mash-paper]
> *Nature Genetics* **51**, 187-195.

## Setup

The simplest and quickest to obtain mashr is to install it from
[CRAN][mashr-cran]:

```R
install.packages("mashr")
```

For Mac and Windows, usually the package is precompiled, so no
additional compiler setup is required. If you have Linux, and in some
other rare cases, you will need to make sure that your R installation
is properly set up to compile packages with C++ source; in particular,
the C++ compiler programs supported by your version of R should be
installed on your computer, and R should be correctly configured to
call these compilers when installing packages from source. For more
information, see the [CRAN documentation][cran-docs].

Alternatively, you may use [devtools][devtools] to install the latest
version of mashr from GitHub:

```R
install.packages("devtools")
devtools::install_github("stephenslab/mashr")
```
   
This command should automatically install any missing dependencies
that are available from CRAN and GitHub. 

## Developer notes

+ When any changes are made to `roxygen2` markup or the C++ code in
the src directory, run `devtools::document()` to update the
[RcppExports.cpp](src/RcppExports.cpp), the package namespaces (see
[NAMESPACE](NAMESPACE)), and the package documentation files (in the
"man" subdirectory),

+ These are the R commands to build the website (make sure you are
connected to Internet while running these commands):

   ```R
   pkgdown::build_site(lazy=TRUE, examples=FALSE)
   ```

+ After editing C++ code in the `src` directory, please use
[uncrustify][uncrustify] to format the code using configuration file
`inst/misc/uncrustify_default.cfg`. For example:

   ```bash
   uncrustify -c uncrustify_default.cfg --replace --no-backup -l CPP mash.cpp
   ```

+ Prior to submitting the package to CRAN, the following modifications
need to be made: (1) remove the `Remotes:` entry in `DESCRIPTION`; (2)
remove the `flash_mash.Rmd` vignette; (3) remove "flashr" from
`Suggests:` in `DESCRIPTION`; (4) Make sure version number is of the
form x.y.z.

[mashr-pkg-for-paper]: http://github.com/stephenslab/mashr-paper
[cran-docs]: https://cran.r-project.org/manuals.html
[mash-paper]: https://doi.org/10.1038/s41588-018-0268-8
[mashr-release]: https://github.com/stephenslab/mashr/releases/tag/v0.2-11
[devtools]: https://github.com/r-lib/devtools
[flashr]: https://github.com/stephenslab/flashr
[vignette-intro]: https://stephenslab.github.io/mashr/articles/intro_mash.html
[vignette-data-driven-cov]: https://stephenslab.github.io/mashr/articles/intro_mash_dd.html
[vignette-non-canonical]: https://stephenslab.github.io/mashr/articles/simulate_noncanon.html
[uncrustify]: http://uncrustify.sourceforge.net
[mashr-cran]: https://cran.r-project.org/package=mashr
