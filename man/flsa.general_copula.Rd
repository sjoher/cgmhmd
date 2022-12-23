\name{flsa.general_copula}
\alias{flsa.general_copula}

\title{flsa.general_copula}

\description{This is a function that minimizes the negative log likelihood over values for the precision matrix \eqn{k \neq 2}.}

\usage{
flsa.general_copula(A,L,lam1,lam2,penalize.diagonal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{List containing \eqn{k} weighted version of \egn{Theta}.}
  \item{L}{number of penalty values to evaluate (1).}
  \item{lam1}{Value for lambda1.}
  \item{lam2}{Value for lambda2.}
  \item{penalize.diagonal}{Choice whether to penalize the diagonal (TRUE) or not (FALSE)}
}

\value{
\item{Z}{Optimized version of the \eqn{k} \egn{Theta}.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2022). Copula graphical models for heterogeneous mixed data. arXiv preprint, arXiv:2210.13140. \cr
}

\author{
  Sjoerd Hermes, Joost van Heerwaarden and Pariya Behrouzi\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
A <- list(matrix(runif(25), 5, 5),matrix(runif(25), 5, 5),matrix(runif(25), 5, 5))
lam1 <- 0.1
lam2 <- 0.5

Z = flsa.general_copula(A,1,lam1,lam2,penalize.diagonal=TRUE)
}
