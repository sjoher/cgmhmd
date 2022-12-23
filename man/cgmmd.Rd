\name{cgmmd}
\alias{cgmmd}

\title{cgmmd}

\description{
This function implements either the Gibbs or approximation method within the Gaussian copula graphical model to estimate the conditional expectation for the data that not follow Gaussianity assumption (e.g. ordinal, discrete, continuous non-Gaussian, or mixed dataset).
}
\usage{
cgmmd(X, method, lambda1, lambda2, ncores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A list containing \egn{k} \eqn{n_k \times p} matrices (\eqn{n_k} is the sample size for group \egn{k} and \eqn{p} is the number of variables}
  \item{method}{Choice between "Gibbs" and "Approximate" indicating which method to use.}
  \item{lambda1}{Vector containing values (in [0,1]) for the sparsity penalization of each \egn{Theta^k}.}
  \item{lambda1}{Vector containing values (in [0,1]) for the similarity penalization between the \egn{Theta^k}.}
  \item{ncores}{Number of cores to be used during parallel computing.}
}

\value{
\item{Z}{ New transformation of the data based on given or default \code{Sigma}.}
\item{ES}{Expectation of covariance matrix( diagonal scaled to 1) of the Gaussian copula graphical model.}
\item{Sigma}{The covariance matrix of the latent variable given the data.}
\item{Theta}{The inverse covariance matrix of the latent variable given the data.}
\item{loglik}{Value of the Log likelihood under the estimated parameters.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2022). Copula graphical models for heterogeneous mixed data. arXiv preprint, arXiv:2210.13140. \cr
}

\author{
  Sjoerd Hermes, Joost van Heerwaarden and Pariya Behrouzi\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
\dontrun{
y <- list(matrix(runif(25), 5, 5),matrix(runif(25), 5, 5),matrix(runif(25), 5, 5))
l1 <- seq(0, 1, 0.1)
l2 <- seq(0, 1, 0.1)
est <- cgmmd(y, "Approximate", l1, l2)
}}

