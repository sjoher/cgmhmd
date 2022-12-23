\name{cutoffs}
\alias{cutoffs}

\title{cutoffs}

\description{
Calculates cut-points of variables with respect to the Gaussian copula.
}

\usage{
cutoffs(y)
}

\arguments{
  \item{y}{
A \eqn{n_k \times p} matrix (\eqn{n_k} is the sample size for group \eqn{k} and \eqn{p} is the number of variables).
}
}

\details{ The relationship between \eqn{j}th variable and \eqn{j}th latent variable is expressed through this set of cut-points.
}

\value{
\item{cutoffs}{ A \eqn{p} by \eqn{(r + 1)} matrix representing the cut-point values under the Gaussian copula, where \eqn{r} defines the number of categories in the dataset.}
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
cutoffs(y[[1]])
}
