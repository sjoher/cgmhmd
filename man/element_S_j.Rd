\name{element_S_j}
\alias{element_S_j}

\title{
element_S_j
}

\description{
Computes the first and second moments of truncated normal random variables}
\usage{
element_S_j(j, lower_upper, mu=0, sigma=1)
}

\arguments{
  \item{j}{Column index (variable \eqn{j} from 1 through \eqn{p})}
  \item{lower_upper}{Matrix containing the lower and upper bands for each data point.}
  \item{mu}{Mean of (standard) normal random variable.}
  \item{sigma}{Standard deviation of (standard) normal random variable.}
}

\value{
\item{EX}{First moment.}
\item{EXX}{Second moment.}
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
lower_upper <- lower.upper(y)
lapply(1:p, function(i){ element_S_j(i, lower_upper))}
}
