\name{marpoly_root_allocations}
\alias{marpoly_root_allocations}
\title{Enumerate MAR root allocations}
\usage{marpoly_root_allocations(ar_coef, tol = 1e-7)}
\arguments{
\item{ar_coef}{Pseudo-causal AR coefficients.}
\item{tol}{Numerical tolerance.}
}
\description{Factorizes a pseudo-causal AR polynomial and allocates inverse roots between causal and noncausal components.}
\value{A data frame with one row per feasible allocation.}
