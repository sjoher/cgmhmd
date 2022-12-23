\name{soft_copula}
\alias{soft_copula}

\title{soft_copula}

\description{
Copmputes soft thresholding function for matrix.}

\usage{
soft_copula(a,lam,penalize.diagonal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Matrix which we want to soft threshold.}
  \item{lam1}{Value for lambda1.}
  \item{penalize.diagonal}{Choice whether to penalize the diagonal (TRUE) or not (FALSE)}
}

\value{
\item{out}{Soft-thresholded version of a.}
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
lam1 <- 0.1

X = soft_copula(y,lam=lam1, penalize.diagonal=TRUE)
}
