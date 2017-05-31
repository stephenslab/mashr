Welcome to Multivariate Adaptive Shrinkage (R package)!

The package implements methods to estimate and test many effects in many conditions (or many effects on many outcomes).

The methods use Empirical Bayes methods to estimate patterns of similarity among conditions, and then exploit
those patterns of similarity among conditions to improve accuracy of effect estimates.
See [Urbut et al](http://biorxiv.org/content/early/2017/05/09/096552) for details of the model and methods.

To install:

```
devtools::install_github("stephenslab/mashr")
library("mashr")
```

See the [Introductory Vignette](vignette/mash_intro.html) for an introduction.

Note that this R package is a refactoring of the code originally used to create results for the paper.
The original package code is [here](http://github.com/stephenslab/mashr-paper).
