## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,comment = "#",fig.width = 5,
                      fig.height = 4,fig.align = "center",
                      eval = FALSE)

## ------------------------------------------------------------------------
#  library(ashr)
#  library(mashr)
#  set.seed(1)
#  simdata = simple_sims(500,5,1)

## ------------------------------------------------------------------------
#  data   = mash_set_data(simdata$Bhat, simdata$Shat)
#  m.1by1 = mash_1by1(data)
#  strong = get_significant_results(m.1by1,0.05)

## ------------------------------------------------------------------------
#  U.pca = cov_pca(data,5,strong)
#  print(names(U.pca))

## ------------------------------------------------------------------------
#  U.ed = cov_ed(data, U.pca, strong)

## ------------------------------------------------------------------------
#  m.ed = mash(data, U.ed)
#  print(get_loglik(m.ed),digits = 10)

## ------------------------------------------------------------------------
#  U.c = cov_canonical(data)
#  m   = mash(data, c(U.c,U.ed))

## ----info----------------------------------------------------------------
#  print(sessionInfo())

