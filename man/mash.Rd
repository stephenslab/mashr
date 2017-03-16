% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mash.R
\name{mash}
\alias{mash}
\title{Apply mash method to data}
\usage{
mash(Bhat, Shat, cov_methods = c("null", "identity", "singletons",
  "equal_effects", "simple_het"), gridmult = sqrt(2), grid = NULL,
  normalizeU = TRUE, prior = "nullbiased", optmethod = c("mixIP"))
}
\arguments{
\item{Bhat}{an n by R matrix of observations (n units in R conditions)}

\item{Shat}{an n by R matrix of standard errors (n units in R conditions)}

\item{cov_methods}{a string indicating what covariance matrix methods to use}

\item{gridmult}{scalar indicating factor by which adjacent grid values should differ; close to 1 for fine grid}

\item{grid}{vector of grid values to use (scaling factors omega in paper)}

\item{prior}{indicates what penalty to use on the likelihood}

\item{optmethod}{name of optimization method to use}
}
\description{
Apply mash method to data
}