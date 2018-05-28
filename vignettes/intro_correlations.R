## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,comment = "#",fig.width = 5,
                      fig.height = 4,fig.align = "center",
                      eval = TRUE)

## ------------------------------------------------------------------------
library(ashr)
library(mashr)
set.seed(1)
simdata = simple_sims(500,5,1)

## ------------------------------------------------------------------------
data   = mash_set_data(simdata$Bhat, simdata$Shat)
V = estimate_null_correlation(data)
data.V = mash_set_data(simdata$Bhat, simdata$Shat, V=V)

## ------------------------------------------------------------------------
U.c = cov_canonical(data.V) 
m.c = mash(data.V, U.c) # fits with correlations because data.V includes correlation information 
print(get_loglik(m.c),digits=10) # log-likelihood of the fit with correlations set to V

## ------------------------------------------------------------------------
m.c.orig = mash(data, U.c) # fits without correlations because data object was set up without correlations
print(get_loglik(m.c.orig),digits=10)

