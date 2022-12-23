\name{lower.upper}
\alias{lower.upper}

\title{
lower.upper
}

\description{
Calculates lower and upper bands for each data point, using a set of cut-points which is obtained from the Gaussian copula.
}

\usage{
lower.upper(y)
}

\arguments{
	\item{y}{ An (\eqn{n_k \times p}) matrix corresponding to the data matrix (\eqn{n_k} is the sample size for group \eqn{k} and \eqn{p} is the number of variables).
}
}

\value{
\item{lower}{
A \eqn{n_k} by \eqn{p} matrix representing the lower band for each data point.
}
\item{upper}{
A \eqn{n_k} by \eqn{p} matrix representing the upper band for each data point.
}
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
lower.upper(y[[1]])
}
