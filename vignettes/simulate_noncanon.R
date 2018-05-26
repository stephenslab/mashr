## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,comment = "#",fig.width = 5,
                      fig.height = 4,fig.align = "center",
                      eval = FALSE)

## ------------------------------------------------------------------------
#  library(ashr)
#  library(mashr)
#  set.seed(1)
#  simdata = simple_sims2(1000,1)
#  true.U1 = cbind(c(1,1,0,0,0),c(1,1,0,0,0),rep(0,5),rep(0,5),rep(0,5))
#  true.U2 = cbind(rep(0,5),rep(0,5),c(0,0,1,1,1),c(0,0,1,1,1),c(0,0,1,1,1))
#  U.true  = list(true.U1 = true.U1, true.U2 = true.U2)

## ---- collapse = TRUE----------------------------------------------------
#  data   = mash_set_data(simdata$Bhat, simdata$Shat)
#  m.1by1 = mash_1by1(data)
#  strong = get_significant_results(m.1by1)
#  U.c    = cov_canonical(data)
#  U.pca  = cov_pca(data,5,strong)
#  U.ed   = cov_ed(data,U.pca,strong)
#  
#  # Computes covariance matrices based on extreme deconvolution,
#  # initialized from PCA.
#  m.c    = mash(data, U.c)
#  m.ed   = mash(data, U.ed)
#  m.c.ed = mash(data, c(U.c,U.ed))
#  m.true = mash(data, U.true)
#  
#  print(get_loglik(m.c),digits = 10)
#  print(get_loglik(m.ed),digits = 10)
#  print(get_loglik(m.c.ed),digits = 10)
#  print(get_loglik(m.true),digits = 10)

