\name{initialize}
\alias{initialize}

\title{initialize}

\description{Initializes parameters to be used in the approximate method algorithm.}

\usage{
initialize(y, ncores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Data.}
  \item{ncores}{Number of cores to be used during parallel computing.}
}
\value{
\item{ES}{Expectation of covariance matrices ( diagonal scaled to 1) of the Gaussian copula graphical model.}
\item{Z}{ New transformation of the data based on given or default \code{Sigma}.}
\item{lower_upper}{ Lower and upper truncation points for the truncated normal distribution.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2022). Copula graphical models for heterogeneous mixed data. arXiv preprint, arXiv:2210.13140. \cr
}

\author{
  Sjoerd Hermes, Joost van Heerwaarden and Pariya Behrouzi\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
y <- list(matrix(runif(25), 5, 5),matrix(runif(25), 5, 5),matrix(runif(25), 5, 5))
ncores <- 7
initialize(y, ncores)
}
