% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NLSS_sum.R
\name{NLSS_sum}
\alias{NLSS_sum}
\title{Summary of the MCMC result for NLSS}
\usage{
NLSS_sum(res, th = 0.95, nstart = 1, nend = 1)
}
\arguments{
\item{res}{result from the function NLSS}
}
\description{
The function summarizes the MCMC result and returns the posterior mean
of A, the posterior mode of S, beta coefficient (frequency of each discrete value of S
among the MCMC samples) and the log-likelihood trace.
}
