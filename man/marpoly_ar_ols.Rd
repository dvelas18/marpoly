\name{marpoly_ar_ols}
\alias{marpoly_ar_ols}
\alias{marpoly_select_p}
\title{Pseudo-causal AR estimation and order selection}
\usage{marpoly_ar_ols(x, p, demean = TRUE, include_mean = FALSE)
marpoly_select_p(x, p_max = 6, criterion = c("BIC", "AIC"), demean = TRUE, include_mean = FALSE)}
\arguments{
\item{x}{Numeric time series.}
\item{p}{AR order.}
\item{demean}{Whether to demean the series.}
\item{include_mean}{Whether to include intercept.}
\item{p_max}{Maximum order.}
\item{criterion}{AIC or BIC.}
}
\description{Fits the pseudo-causal AR(p) representation by OLS and selects p by information criteria.}
\value{A \code{marpoly_ar} or \code{marpoly_ar_selection} object.}
