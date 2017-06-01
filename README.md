# mashr: Multivariate Adaptive Shrinkage in R

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
and [Simulation with non-canonical matrices](simulate_noncanon.html).

## Setup

These are the necessary steps to install the most recent version of
the mashr package:

1. In R, install these three R packages from CRAN:

   ```R
   install.packages(c("assertthat","mvtnorm","rmeta"))
   ```

2. Optionally, install MOSEK and the Rmosek package, for faster
   optimization in the `ashr` package. See the
   [ashr Github repository](https://github.com/stephens999/ashr) for
   details.

3. Install the [ExtremeDeconvolution R package](https://github.com/jobovy/extreme-deconvolution#installation). Note that you will need to link to the
   [GNU Scientific Library](https://www.gnu.org/software/gsl) to
   build this package.

4. Once you have installed all these packages, you can install and
   load the most recent version of `mashr` available on Github. This
   command will automatically retrieve and install version 2.1-19 of
   the `ashr` package released on Github.

   ```R
   library(devtools)
   install_github("stephenslab/mashr")
   library(mashr)
   ```
