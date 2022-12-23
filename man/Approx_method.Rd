\name{Approx_method}
\alias{Approx_method}

\title{Approx_method}

\description{
This function implements the approximation method within the Gaussian copula graphical model to estimate the conditional expectation for the data that not follow Gaussianity assumption (e.g. ordinal, discrete, continuous non-Gaussian, or mixed dataset).
}

\usage{
Approx_method(y, Z, ES=NULL, lambdas, lower_upper=NULL, combination = 1, em_tol=.001, em_iter=10, ncores)
}

\arguments{
	\item{y}{A list containing \egn{k} \eqn{n_k \times p} matrices (\eqn{n_k} is the sample size for group \egn{k} and \eqn{p} is the number of variables).
	\item{Z}{ A list containing \egn{k} \eqn{n_k \times p} matrices which are transformations of the data via the Gaussian copula. If \code{Z = NULL} internally calculates an initial value for \eqn{Z}.}
	\item{ES}{Expectation of covariance matrices ( diagonal scaled to 1) of the Gaussian copula graphical model.}
	\item{lambdas}{A matrix containing values of two (non-negative) regularization parameters to calculate \code{Sigma^k}. lambdas=0 means no regularization.}
	\item{lower_upper}{Two matrices containing respectively the lower and upper bounds for each data point.}
	\item{combination}{Number from 1:nrow(lambdas) over which the parameters are computed.}
	\item{em_tol}{Minimum tolerance value at which the EM algorithm terminates.}
	\item{em_iter}{Maximum number of iterations used by the EM algorithm.}
	\item{ncores}{If \code{ncores = NULL}, the algorithm internally detects number of available cores and run the calculations in parallel on (available cores - 1). Typical usage is to fix \code{ncores = 1} when \eqn{p} is small \eqn{( p < 500 )}, and \code{ncores = NULL} when \eqn{p} is large.}
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
ncores <- 7
ini <- initialize(y, ncores = ncores)
Z	= ini$Z
ES	= ini$ES
lower_upper = ini$lower_upper
lambda1 <- seq(0,1,0.1)
lambda2 <- seq(0,1,0.1)
Approx_method(y, Z, ES, lambdas = c(lambda1, lambda2), lower_upper, combination = 1, em_tol=.001, em_iter=10, ncores)
}
}
