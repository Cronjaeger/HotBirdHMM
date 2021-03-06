% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HMM.R
\name{sample_XandY}
\alias{sample_XandY}
\title{Sample full model}
\usage{
sample_XandY(rho, theta, M, L, sample_x0 = NULL)
}
\arguments{
\item{rho, }{rathe of recombination between ajacent sites; affects transition
rate of hidden states.}

\item{theta}{Rate of mutation; affects distribution of Y|X.}

\item{M}{number of hidden states of X.}

\item{L}{desired length of sampled sequences.}

\item{sample_x0}{a function which generates samples the initial state of the
hidden variables X.}
}
\value{
a list with names “X” and “Y”
X is a vector of length \code{L} encoding the hidden states
Y is a vector of length \code{L} encoding the emissions.
}
\description{
Samples both hidden states and emissions
}
