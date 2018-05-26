## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,comment = "#",fig.width = 5,
                      fig.height = 4,fig.align = "center")

## ------------------------------------------------------------------------
library(ashr)
library(mashr)
set.seed(1)
simdata = simple_sims(500,5,1)

## ------------------------------------------------------------------------
data = mash_set_data(simdata$Bhat, simdata$Shat)

## ------------------------------------------------------------------------
U.c = cov_canonical(data)  
print(names(U.c))

## ------------------------------------------------------------------------
m.c = mash(data, U.c)

## ------------------------------------------------------------------------
head(get_lfsr(m.c))
head(get_pm(m.c))
head(get_psd(m.c))

## ------------------------------------------------------------------------
head(get_significant_results(m.c))
print(length(get_significant_results(m.c)))

## ------------------------------------------------------------------------
print(head(get_significant_results(m.c, conditions=1)))

## ------------------------------------------------------------------------
print(get_pairwise_sharing(m.c)) 

## ------------------------------------------------------------------------
print(get_pairwise_sharing(m.c, factor=0))

## ------------------------------------------------------------------------
print(get_loglik(m.c))

## ------------------------------------------------------------------------
print(get_estimated_pi(m.c))
barplot(get_estimated_pi(m.c),las = 2)

## ------------------------------------------------------------------------
mash_plot_meta(m.c,get_significant_results(m.c)[1])

## ----info----------------------------------------------------------------
print(sessionInfo())

