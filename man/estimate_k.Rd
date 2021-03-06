% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_functions.R
\name{estimate_k}
\alias{estimate_k}
\title{estimate_k: estimate K using the eigengap or multimodality gap heuristics}
\usage{
estimate_k(A2, maxk = 10, showplots = TRUE)
}
\arguments{
\item{A2}{Data frame or matrix: a similarity matrix}

\item{maxk}{Numerical value: maximum number of K to be considered}

\item{showplots}{Character value: whether to show the plot on the screen}
}
\value{
A data frame containing the eigenvalues and dip-test statistics of the 
eigenvectors of the graph Laplacian
}
\description{
This function will try to estimate K given a similarity matrix. Generally the
maximum eigengap is preferred, but on some data examining the distribution
of the eigenvectors as in the multimodality gap heuristic may be beneficial.
}
\examples{
k_test <- estimate_k(missl[[1]])
}
